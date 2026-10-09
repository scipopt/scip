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

/**@file   setters.c
 * @brief  unit test for checking setters
 * @author Franziska Schloesser
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"

/* UNIT TEST */

#include "include/scip_test.h"

/* GLOBAL VARIABLES */
static SCIP* scip;

/** create unbounded problem */
static
void initProb(void)
{
   SCIP_CONS* cons;
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 1.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", -SCIPinfinity(scip), SCIPinfinity(scip), -1.0, SCIP_VARTYPE_INTEGER) );

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
}

static
void setAndTestObjsense(SCIP_OBJSENSE sense)
{
   SCIP_CALL( SCIPsetObjsense(scip, sense) );

   TEST_ASSERT_EQUAL(sense, SCIPgetObjsense(scip));
}

/* TEST SUITES */
/** setup of test run */
static
void setup(void)
{
   /* initialize SCIP */
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );
   TEST_ASSERT_NOT_NULL(scip);

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

/** deinitialization method */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

void setUp(void)
{
   setup();
}

void tearDown(void)
{
   teardown();
}

/* TESTS */
void test_setters_setProbNameTest(void)
{
   char name[SCIP_MAXSTRLEN] = "new_na-me";
   SCIP_CALL( SCIPsetProbName(scip, name) );

   TEST_ASSERT_EQUAL_STRING(name, SCIPgetProbName(scip));
}

void test_setters_setObjsenseTest(void)
{
   setAndTestObjsense(SCIP_OBJSENSE_MAXIMIZE);
   setAndTestObjsense(SCIP_OBJSENSE_MINIMIZE);

   /* add some vars and test again */
   initProb();

   setAndTestObjsense(SCIP_OBJSENSE_MAXIMIZE);
   setAndTestObjsense(SCIP_OBJSENSE_MINIMIZE);
}

/** macro to check the set methods for priorities
 *
 *
 */
#define TEST_PRIORITY(pointertype, getobjs, getnobjs, setpriority, getpriority)           \
   do                                                                                     \
   {                                                                                      \
      int i;                                                                              \
      int priority;                                                                       \
      pointertype objs;                                                                   \
      int nobjs;                                                                          \
                                                                                          \
      objs = getobjs(scip);                                                               \
      nobjs = getnobjs(scip);                                                             \
                                                                                          \
      /* set and test priorities */                                                       \
      priority = 11;                                                                      \
      for( i = 0; i < nobjs; ++i )                                                        \
      {                                                                                   \
         SCIP_CALL(setpriority(scip, objs[i], priority));                                 \
         TEST_ASSERT_EQUAL(priority, getpriority(objs[i]));                                    \
      }                                                                                   \
                                                                                          \
      /* add some vars and test priorities again */                                       \
      initProb();                                                          \
                                                                                          \
      for( i=0; i < nobjs; ++i)                                                           \
      {                                                                                   \
         TEST_ASSERT_EQUAL(priority, getpriority(objs[i]));                                    \
      }                                                                                   \
                                                                                          \
                                                                                          \
      /* set and test priorities */                                                       \
      priority = 12;                                                                      \
      for( i = 0; i < nobjs; ++i )                                                        \
      {                                                                                   \
         SCIP_CALL(setpriority(scip, objs[i], priority));                                 \
         TEST_ASSERT_EQUAL(priority, getpriority(objs[i]));                                    \
      }                                                                                   \
   } while(FALSE)                                                                         \

void test_setters_setConflicthdlrPriorityTest(void)
{
   TEST_PRIORITY( SCIP_CONFLICTHDLR**, SCIPgetConflicthdlrs, SCIPgetNConflicthdlrs, SCIPsetConflicthdlrPriority, SCIPconflicthdlrGetPriority );
}

void test_setters_setPresolPriorityTest(void)
{
   TEST_PRIORITY( SCIP_PRESOL**, SCIPgetPresols, SCIPgetNPresols, SCIPsetPresolPriority, SCIPpresolGetPriority );
}

void test_setters_setPricerPriorityTest(void)
{
   TEST_PRIORITY( SCIP_PRICER**, SCIPgetPricers, SCIPgetNPricers, SCIPsetPricerPriority, SCIPpricerGetPriority );
}

void test_setters_setRelaxPriorityTest(void)
{
   TEST_PRIORITY( SCIP_RELAX**, SCIPgetRelaxs, SCIPgetNRelaxs, SCIPsetRelaxPriority, SCIPrelaxGetPriority );
}

void test_setters_setSepaPriorityTest(void)
{
   TEST_PRIORITY( SCIP_SEPA**, SCIPgetSepas, SCIPgetNSepas, SCIPsetSepaPriority, SCIPsepaGetPriority );
}

void test_setters_setPropPriorityTest(void)
{
   TEST_PRIORITY( SCIP_PROP**, SCIPgetProps, SCIPgetNProps, SCIPsetPropPriority, SCIPpropGetPriority );
}

void test_setters_setHeurPriorityTest(void)
{
   TEST_PRIORITY( SCIP_HEUR**, SCIPgetHeurs, SCIPgetNHeurs, SCIPsetHeurPriority, SCIPheurGetPriority );
}

void test_setters_setNodeselStdPriorityTest(void)
{
   TEST_PRIORITY( SCIP_NODESEL**, SCIPgetNodesels, SCIPgetNNodesels, SCIPsetNodeselStdPriority, SCIPnodeselGetStdPriority );
}

void test_setters_setNodeselMemsavePriorityTest(void)
{
   TEST_PRIORITY( SCIP_NODESEL**, SCIPgetNodesels, SCIPgetNNodesels, SCIPsetNodeselMemsavePriority, SCIPnodeselGetMemsavePriority );
}

void test_setters_setBranchrulePriorityTest(void)
{
   TEST_PRIORITY( SCIP_BRANCHRULE**, SCIPgetBranchrules, SCIPgetNBranchrules, SCIPsetBranchrulePriority, SCIPbranchruleGetPriority );
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_setters_setProbNameTest);
   RUN_TEST(test_setters_setObjsenseTest);
   RUN_TEST(test_setters_setConflicthdlrPriorityTest);
   RUN_TEST(test_setters_setPresolPriorityTest);
   RUN_TEST(test_setters_setPricerPriorityTest);
   RUN_TEST(test_setters_setRelaxPriorityTest);
   RUN_TEST(test_setters_setSepaPriorityTest);
   RUN_TEST(test_setters_setPropPriorityTest);
   RUN_TEST(test_setters_setHeurPriorityTest);
   RUN_TEST(test_setters_setNodeselStdPriorityTest);
   RUN_TEST(test_setters_setNodeselMemsavePriorityTest);
   RUN_TEST(test_setters_setBranchrulePriorityTest);
   return UNITY_END();
}
