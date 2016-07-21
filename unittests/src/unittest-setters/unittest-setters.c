/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  unit test for checking setters on scip.c
 *
 * DONE:
 *
 * SCIP_RETCODE SCIPsetProbName
 * SCIP_RETCODE SCIPsetObjsense
 * SCIP_RETCODE SCIPsetConflicthdlrPriority
 * SCIP_RETCODE SCIPsetPresolPriority
 * SCIP_RETCODE SCIPsetPricerPriority  @todo add a pricer
 * SCIP_RETCODE SCIPsetRelaxPriority  @todo add a relax
 * SCIP_RETCODE SCIPsetPropPriority
 * SCIP_RETCODE SCIPsetHeurPriority
 * SCIP_RETCODE SCIPsetNodeselStdPriority
 * SCIP_RETCODE SCIPsetNodeselMemsavePriority
 * SCIP_RETCODE SCIPsetBranchrulePriority
 *
 * @todo
 * SCIP_RETCODE SCIPsetMessagehdlr
 * void SCIPsetMessagehdlrLogfile(
 * void SCIPsetMessagehdlrQuiet(
 * SCIP_RETCODE SCIPsetParam(
 * SCIP_RETCODE SCIPsetBoolParam(
 * SCIP_RETCODE SCIPsetIntParam(
 * SCIP_RETCODE SCIPsetLongintParam(
 * SCIP_RETCODE SCIPsetRealParam(
 * SCIP_RETCODE SCIPsetCharParam(
 * SCIP_RETCODE SCIPsetStringParam(
 * SCIP_RETCODE SCIPsetEmphasis
 * SCIP_RETCODE SCIPsetSubscipsOff
 * SCIP_RETCODE SCIPsetHeuristics
 * SCIP_RETCODE SCIPsetPresolving
 * SCIP_RETCODE SCIPsetSeparating
 * SCIP_RETCODE SCIPsetReaderCopy
 * SCIP_RETCODE SCIPsetReaderFree
 * SCIP_RETCODE SCIPsetReaderRead
 * SCIP_RETCODE SCIPsetReaderWrite
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include "scip/scipdefplugins.h"
#include <string.h>

/** macro to check the return of tests
 *
 *  @note assumes the existence of SCIP_RETCODE retcode
 */
#define CHECK_TEST(x)                            \
   do                                            \
   {                                             \
      retcode = (x);                             \
      if( retcode != SCIP_OKAY )                 \
      {                                          \
         printf("Unit test " #x " failed\n");    \
         SCIPprintError(retcode);                \
         return -1;                              \
      }                                          \
   }                                             \
   while( FALSE )

/** macro to check the set methods for priorities
 *
 *
 */
#define TEST_PRIORITY(pointertype, getobjs, getnobjs, setpriority, getpriority)           \
   do                                                                                     \
   {                                                                                      \
      int i;                                                                              \
      int priority;                                                                       \
      SCIP* scip;                                                                         \
      pointertype objs;                                                                   \
      int nobjs;                                                                          \
                                                                                          \
      SCIP_CALL( SCIPcreate(&scip) );                                                     \
                                                                                          \
      /* include default plugins */                                                       \
      SCIP_CALL( SCIPincludeDefaultPlugins(scip) );                                       \
                                                                                          \
      /* create problem */                                                                \
      SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );                                  \
                                                                                          \
      objs = getobjs(scip);                                                               \
      nobjs = getnobjs(scip);                                                             \
                                                                                          \
      /* set and test priorities */                                                       \
      priority = 11;                                                                      \
      for( i = 0; i < nobjs; ++i )                                                        \
      {                                                                                   \
         SCIP_CALL( setpriority(scip, objs[i], priority) );                               \
         if( priority != getpriority(objs[i]) )                                           \
         {                                                                                \
            return SCIP_ERROR;                                                            \
         }                                                                                \
      }                                                                                   \
                                                                                          \
      /* add some vars and test priorities again */                                       \
      SCIP_CALL( initProb(scip) );                                                        \
                                                                                          \
      for( i=0; i < nobjs; ++i)                                                           \
      {                                                                                   \
         if( priority != getpriority(objs[i]) )                                           \
         {                                                                                \
            return SCIP_ERROR;                                                            \
         }                                                                                \
      }                                                                                   \
                                                                                          \
                                                                                          \
      /* set and test priorities */                                                       \
      priority = 12;                                                                      \
      for( i = 0; i < nobjs; ++i )                                                        \
      {                                                                                   \
         SCIP_CALL( setpriority(scip, objs[i], priority) );                               \
         if( priority != getpriority(objs[i]) )                                           \
         {                                                                                \
            return SCIP_ERROR;                                                            \
         }                                                                                \
      }                                                                                   \
   }                                                                                      \
   while(FALSE)

/** macro to set and check parameters on the constraints
 * @note assumes the existence of SCIP_CONS* cons and SCIP* scip
 */
#define CONS_TEST(consset, iscons)        \
   do{                                    \
      conset(scip, cons, true);           \
      if( iscons(scip, cons) == FALSE )   \
         return SCIP_ERROR;               \
                                          \
      conset(scip, cons, false);          \
      if( iscons(scip,cons) == TRUE )     \
         return SCIP_ERROR;               \
   }                                      \
   while(FALSE)

/*
 * Local methods for tests
 */

/** create unbounded problem */
static
SCIP_RETCODE initProb(
   SCIP*                 scip                /**< SCIP instance */
   )
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

   return SCIP_OKAY;
}

static
SCIP_RETCODE setAndTestObjsense(SCIP* scip, SCIP_OBJSENSE sense)
{
   SCIP_CALL( SCIPsetObjsense(scip, sense) );

   if( sense != SCIPgetObjsense(scip) )
   {
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/*
 * TESTS
 */

/** test setProbName */
static
SCIP_RETCODE setProbNameTest(void)
{
   SCIP* scip;
   char name[SCIP_MAXSTRLEN];

   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "originalName") );

   /* change name and test */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "new_na-me");
   SCIP_CALL( SCIPsetProbName(scip, name) );

   if( strcmp(name, SCIPgetProbName(scip)) )
   {
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** test setObjsense */
static
SCIP_RETCODE setObjsenseTest(void)
{
   SCIP* scip;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   SCIP_CALL( setAndTestObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   SCIP_CALL( setAndTestObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* add some vars and test again */
   SCIP_CALL( initProb(scip) );

   SCIP_CALL( setAndTestObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   SCIP_CALL( setAndTestObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   return SCIP_OKAY;
}

/** test setConflicthdlrPriority */
static
SCIP_RETCODE setConflicthdlrPriorityTest(void)
{
   TEST_PRIORITY( SCIP_CONFLICTHDLR**, SCIPgetConflicthdlrs, SCIPgetNConflicthdlrs, SCIPsetConflicthdlrPriority, SCIPconflicthdlrGetPriority );

   return SCIP_OKAY;
}


/** test setPresolPriority */
static
SCIP_RETCODE setPresolPriorityTest(void)
{
   TEST_PRIORITY( SCIP_PRESOL**, SCIPgetPresols, SCIPgetNPresols, SCIPsetPresolPriority, SCIPpresolGetPriority );

   return SCIP_OKAY;
}

/** test setPricerPriority */
static
SCIP_RETCODE setPricerPriorityTest(void)
{
   TEST_PRIORITY( SCIP_PRICER**, SCIPgetPricers, SCIPgetNPricers, SCIPsetPricerPriority, SCIPpricerGetPriority );

   return SCIP_OKAY;
}

/** test setRelaxPriority */
static
SCIP_RETCODE setRelaxPriorityTest(void)
{
   TEST_PRIORITY( SCIP_RELAX**, SCIPgetRelaxs, SCIPgetNRelaxs, SCIPsetRelaxPriority, SCIPrelaxGetPriority );

   return SCIP_OKAY;
}


/** test setSepaPriority */
static
SCIP_RETCODE setSepaPriorityTest(void)
{
   TEST_PRIORITY( SCIP_SEPA**, SCIPgetSepas, SCIPgetNSepas, SCIPsetSepaPriority, SCIPsepaGetPriority );

   return SCIP_OKAY;
}

/** test setPropPriority */
static
SCIP_RETCODE setPropPriorityTest(void)
{
   TEST_PRIORITY( SCIP_PROP**, SCIPgetProps, SCIPgetNProps, SCIPsetPropPriority, SCIPpropGetPriority );

   return SCIP_OKAY;
}


/** test setHeurPriority */
static
SCIP_RETCODE setHeurPriorityTest(void)
{
   TEST_PRIORITY( SCIP_HEUR**, SCIPgetHeurs, SCIPgetNHeurs, SCIPsetHeurPriority, SCIPheurGetPriority );

   return SCIP_OKAY;
}



/** test setNodeselStdPriority */
static
SCIP_RETCODE setNodeselStdPriorityTest(void)
{
   TEST_PRIORITY( SCIP_NODESEL**, SCIPgetNodesels, SCIPgetNNodesels, SCIPsetNodeselStdPriority, SCIPnodeselGetStdPriority );

   return SCIP_OKAY;
}

/** test setNodeselMemsavePriority */
static
SCIP_RETCODE setNodeselMemsavePriorityTest(void)
{
   TEST_PRIORITY( SCIP_NODESEL**, SCIPgetNodesels, SCIPgetNNodesels, SCIPsetNodeselMemsavePriority, SCIPnodeselGetMemsavePriority );

   return SCIP_OKAY;
}

/** test setNodeselMemsavePriority */
static
SCIP_RETCODE setBranchrulePriorityTest(void)
{
   TEST_PRIORITY( SCIP_BRANCHRULE**, SCIPgetBranchrules, SCIPgetNBranchrules, SCIPsetBranchrulePriority, SCIPbranchruleGetPriority );

   return SCIP_OKAY;
}

/** test setConsInitial */

/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-setters ===========\n");
   printf("=opt=  unittest-setters 0\n\n");

   CHECK_TEST( setProbNameTest() );
   CHECK_TEST( setObjsenseTest() );
   CHECK_TEST( setConflicthdlrPriorityTest() );
   CHECK_TEST( setPresolPriorityTest() );
   CHECK_TEST( setPricerPriorityTest() );
   CHECK_TEST( setRelaxPriorityTest() );
   CHECK_TEST( setSepaPriorityTest() );
   CHECK_TEST( setPropPriorityTest() );
   CHECK_TEST( setHeurPriorityTest()  );
   CHECK_TEST( setNodeselStdPriorityTest() );
   CHECK_TEST( setNodeselMemsavePriorityTest() );
   CHECK_TEST( setBranchrulePriorityTest() );

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return 0;
}
