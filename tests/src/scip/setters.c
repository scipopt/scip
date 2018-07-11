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

   cr_assert_eq(sense, SCIPgetObjsense(scip));
}

/* TEST SUITES */
/** setup of test run */
static
void setup(void)
{
   /* initialize SCIP */
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );
   cr_assert_not_null(scip);

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

/** deinitialization method */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

TestSuite(setters, .init = setup, .fini = teardown);

/* TESTS */
Test(setters, setProbNameTest)
{
   char name[SCIP_MAXSTRLEN] = "new_na-me";
   SCIP_CALL( SCIPsetProbName(scip, name) );

   cr_assert_str_eq(name, SCIPgetProbName(scip));
}

Test(setters, setObjsenseTest)
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
         cr_assert_eq(priority, getpriority(objs[i]));                                    \
      }                                                                                   \
                                                                                          \
      /* add some vars and test priorities again */                                       \
      initProb();                                                          \
                                                                                          \
      for( i=0; i < nobjs; ++i)                                                           \
      {                                                                                   \
         cr_assert_eq(priority, getpriority(objs[i]));                                    \
      }                                                                                   \
                                                                                          \
                                                                                          \
      /* set and test priorities */                                                       \
      priority = 12;                                                                      \
      for( i = 0; i < nobjs; ++i )                                                        \
      {                                                                                   \
         SCIP_CALL(setpriority(scip, objs[i], priority));                                 \
         cr_assert_eq(priority, getpriority(objs[i]));                                    \
      }                                                                                   \
   } while(FALSE)                                                                         \

Test(setters, setConflicthdlrPriorityTest)
{
   TEST_PRIORITY( SCIP_CONFLICTHDLR**, SCIPgetConflicthdlrs, SCIPgetNConflicthdlrs, SCIPsetConflicthdlrPriority, SCIPconflicthdlrGetPriority );
}

Test(setters, setPresolPriorityTest)
{
   TEST_PRIORITY( SCIP_PRESOL**, SCIPgetPresols, SCIPgetNPresols, SCIPsetPresolPriority, SCIPpresolGetPriority );
}

Test(setters, setPricerPriorityTest)
{
   TEST_PRIORITY( SCIP_PRICER**, SCIPgetPricers, SCIPgetNPricers, SCIPsetPricerPriority, SCIPpricerGetPriority );
}

Test(setters, setRelaxPriorityTest)
{
   TEST_PRIORITY( SCIP_RELAX**, SCIPgetRelaxs, SCIPgetNRelaxs, SCIPsetRelaxPriority, SCIPrelaxGetPriority );
}

Test(setters, setSepaPriorityTest)
{
   TEST_PRIORITY( SCIP_SEPA**, SCIPgetSepas, SCIPgetNSepas, SCIPsetSepaPriority, SCIPsepaGetPriority );
}

Test(setters, setPropPriorityTest)
{
   TEST_PRIORITY( SCIP_PROP**, SCIPgetProps, SCIPgetNProps, SCIPsetPropPriority, SCIPpropGetPriority );
}

Test(setters, setHeurPriorityTest)
{
   TEST_PRIORITY( SCIP_HEUR**, SCIPgetHeurs, SCIPgetNHeurs, SCIPsetHeurPriority, SCIPheurGetPriority );
}

Test(setters, setNodeselStdPriorityTest)
{
   TEST_PRIORITY( SCIP_NODESEL**, SCIPgetNodesels, SCIPgetNNodesels, SCIPsetNodeselStdPriority, SCIPnodeselGetStdPriority );
}

Test(setters, setNodeselMemsavePriorityTest)
{
   TEST_PRIORITY( SCIP_NODESEL**, SCIPgetNodesels, SCIPgetNNodesels, SCIPsetNodeselMemsavePriority, SCIPnodeselGetMemsavePriority );
}

Test(setters, setBranchrulePriorityTest)
{
   TEST_PRIORITY( SCIP_BRANCHRULE**, SCIPgetBranchrules, SCIPgetNBranchrules, SCIPsetBranchrulePriority, SCIPbranchruleGetPriority );
}
