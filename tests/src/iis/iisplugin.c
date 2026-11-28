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

/**@file   iisplugin.c
 * @brief  unit test for checking iis functionality of scip.c
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scip_iisfinder.h"
#include "scip/scip_prob.h"
#include "scip/scipdefplugins.h"

#include "include/scip_test.h"
#include <string.h>

/** GLOBAL VARIABLES **/
static SCIP* scip;

/** TEST SUITES **/
static
void setup(void)
{
   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}


/** reads an instance that is located in the directory of this file */
static
void readTestProb(
   const char*           probname            /**< file name of the instance */
   )
{
   char filename[SCIP_MAXSTRLEN];

   TESTsetTestfilename(filename, __FILE__, probname);
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );
}


/* test that the IIS functionality works */
void test_iisplugin_valid(void)
{
   SCIP_IIS* iis;

   readTestProb("test_infeasible.lp");
   SCIP_CALL( SCIPsolve(scip) );
   iis = SCIPgetIIS(scip);
   /** ensure that the original problem is infeasible */
   SOFT_ASSERT_EQUAL(SCIPgetStatus(scip), SCIP_STATUS_INFEASIBLE, "got status %d, expected %d", SCIPgetStatus(scip), SCIP_STATUS_INFEASIBLE);
   /** ensure that the iis does not yet exist and is therefore invalid */
   SOFT_ASSERT_EQUAL(SCIPiisIsSubscipInfeasible(iis), FALSE, "iis is valid before doing any computations");
   SCIP_CALL( SCIPgenerateIIS(scip) );
   /** ensure that the iis exists and is therefore valid */
   SOFT_ASSERT_EQUAL(SCIPiisIsSubscipInfeasible(iis), TRUE, "iis is not valid");

   SCIP* subscip = SCIPiisGetSubscip(iis);
   int nOrigVars = SCIPgetNOrigVars(scip);
   SOFT_ASSERT(nOrigVars > SCIPgetNOrigConss(subscip), "original has more variables than iis");

   /** test a removed and a preserved variable */
   SCIP_VAR** origVars = SCIPgetOrigVars(scip);
   SCIP_VAR* x2 = NULL;
   SCIP_VAR* x6 = NULL;
   for( int i = 0; i < nOrigVars; ++i )
   {
      SCIP_VAR* var = origVars[i];
      if( strcmp(SCIPvarGetName(var), "x2") == 0 )
         x2 = var;
      else if( strcmp(SCIPvarGetName(var), "x6") == 0 )
         x6 = var;
   }
   TEST_ASSERT(x2 != NULL);
   TEST_ASSERT(x6 != NULL);
   SCIP_VAR* x2IIS = SCIPiisGetSubscipVar(iis, x2);
   TEST_ASSERT(x2IIS != NULL);
   TEST_ASSERT(strcmp(SCIPvarGetName(x2), SCIPvarGetName(x2IIS)) == 0);
   SCIP_VAR* x6IIS = SCIPiisGetSubscipVar(iis, x6);
   TEST_ASSERT(x6IIS == NULL);

   /** test a removed and preserved constraint */
   int nOrigConss = SCIPgetNOrigConss(scip);
   SCIP_CONS** origConss = SCIPgetOrigConss(scip);
   SCIP_CONS* c1 = NULL;
   SCIP_CONS* c2 = NULL;
   SCIP_CONS* c5 = NULL;
   for( int i = 0; i < nOrigConss; ++i )
   {
      SCIP_CONS* cons = origConss[i];
      if( strcmp(SCIPconsGetName(cons), "c1") == 0 )
         c1 = cons;
      else if( strcmp(SCIPconsGetName(cons), "c2") == 0 )
         c2 = cons;
      else if( strcmp(SCIPconsGetName(cons), "c5") == 0 )
         c5 = cons;
   }
   TEST_ASSERT(c1 != NULL);
   TEST_ASSERT(c2 != NULL);
   TEST_ASSERT(c5 != NULL);
   SCIP_CONS* c1IIS = SCIPiisGetSubscipCons(iis, c1);
   SCIP_CONS* c2IIS = SCIPiisGetSubscipCons(iis, c2);
   SCIP_CONS* c5IIS = SCIPiisGetSubscipCons(iis, c5);
   TEST_ASSERT(c1IIS == NULL);
   TEST_ASSERT(c2IIS != NULL || c5IIS != NULL);
   TEST_ASSERT(c2IIS == NULL || strcmp(SCIPconsGetName(c2), SCIPconsGetName(c2IIS)) == 0);
   TEST_ASSERT(c5IIS == NULL || strcmp(SCIPconsGetName(c5), SCIPconsGetName(c5IIS)) == 0);
}

/* test that the linear constraint of a deleted indicator constraint does not survive in the IIS */
void test_iisplugin_indicator(void)
{
   SCIP_IIS* iis;

   readTestProb("test_indicator.cip");
   SCIP_CALL( SCIPsolve(scip) );
   /** ensure that the original problem is infeasible */
   SOFT_ASSERT_EQUAL(SCIPgetStatus(scip), SCIP_STATUS_INFEASIBLE, "got status %d, expected %d", SCIPgetStatus(scip), SCIP_STATUS_INFEASIBLE);
   SCIP_CALL( SCIPsetBoolParam(scip, "iis/greedy/additive", FALSE) );
   SCIP_CALL( SCIPsetIntParam(scip, "iis/greedy/maxbatchsize", 1) );
   SCIP_CALL( SCIPgenerateIIS(scip) );
   iis = SCIPgetIIS(scip);
   SOFT_ASSERT_EQUAL(SCIPiisIsSubscipInfeasible(iis), TRUE, "iis is not valid");
   SOFT_ASSERT_EQUAL(SCIPiisIsSubscipIrreducible(iis), TRUE, "iis is not irreducible");

   /** collect the constraints of the instance */
   int nOrigConss = SCIPgetNOrigConss(scip);
   SCIP_CONS** origConss = SCIPgetOrigConss(scip);
   SCIP_CONS* linear1 = NULL;
   SCIP_CONS* linear2 = NULL;
   SCIP_CONS* linear3 = NULL;
   SCIP_CONS* indicator1 = NULL;
   for( int i = 0; i < nOrigConss; ++i )
   {
      SCIP_CONS* cons = origConss[i];
      if( strcmp(SCIPconsGetName(cons), "linear1") == 0 )
         linear1 = cons;
      else if( strcmp(SCIPconsGetName(cons), "linear2") == 0 )
         linear2 = cons;
      else if( strcmp(SCIPconsGetName(cons), "linear3") == 0 )
         linear3 = cons;
      else if( strcmp(SCIPconsGetName(cons), "indicator1") == 0 )
         indicator1 = cons;
   }
   TEST_ASSERT(linear1 != NULL);
   TEST_ASSERT(linear2 != NULL);
   TEST_ASSERT(linear3 != NULL);
   TEST_ASSERT(indicator1 != NULL);

   /** the only IIS consists of linear2 and linear3 */
   SOFT_ASSERT(SCIPiisGetSubscipCons(iis, linear1) == NULL, "linear1 is still part of the iis");
   SOFT_ASSERT(SCIPiisGetSubscipCons(iis, linear2) != NULL, "linear2 is not part of the iis");
   SOFT_ASSERT(SCIPiisGetSubscipCons(iis, linear3) != NULL, "linear3 is not part of the iis");
   SOFT_ASSERT(SCIPiisGetSubscipCons(iis, indicator1) == NULL, "indicator1 is still part of the iis");
   SOFT_ASSERT_EQUAL(SCIPgetNOrigConss(SCIPiisGetSubscip(iis)), 2, "iis has %d constraints, expected 2", SCIPgetNOrigConss(SCIPiisGetSubscip(iis)));
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_iisplugin_valid);
   RUN_TEST(test_iisplugin_indicator);
   return UNITY_END();
}
