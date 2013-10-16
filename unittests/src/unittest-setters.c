/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
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
 *
 * EASY
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

/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   CHECK_TEST( setProbNameTest() );
   CHECK_TEST( setObjsenseTest() );

   printf("All tests passed\n");
   return 0;
}
