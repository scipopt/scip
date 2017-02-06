/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   unittest-presol.c
 * @brief  unit test for checking setters on scip.c
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include "presol_unittest.h"
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



/* local methods */


/** create bounded problem */
static
SCIP_RETCODE initProb(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_CONS* cons;
   SCIP_Real vals[2];
   SCIP_VAR* vars[2];
   char name[SCIP_MAXSTRLEN];

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0, 3, -1.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0, 3, -1.0, SCIP_VARTYPE_INTEGER) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );

   /* add a constraint x + y <= 2 */
   vals[0] = 1.0;
   vals[1] = 1.0;
   vars[0] = xvar;
   vars[1] = yvar;
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "1. Cons");
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, 0.0, 2.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );

   return SCIP_OKAY;
}


/* Check methods */

/*
DONE:
SCIPpresolGetName
SCIPpresolGetDesc
SCIPpresolGetPriority
SCIPpresolIsInitialized

TODO:
SCIPpresolGetData
SCIPpresolWasDelayed
SCIPpresolGetSetupTime
SCIPpresolGetTime
SCIPpresolGetNFixedVars
SCIPpresolGetNAggrVars
SCIPpresolGetNChgVarTypes
SCIPpresolGetNChgBds
SCIPpresolGetNAddHoles
SCIPpresolGetNDelConss
SCIPpresolGetNAddConss
SCIPpresolGetNUpgdConss
SCIPpresolGetNChgCoefs
SCIPpresolGetNChgSides
SCIPpresolGetNCalls
SCIPpresolGetMaxrounds
SCIPpresolGetTiming
SCIPpresolSetTiming
*/



static
SCIP_RETCODE checkPresolGetName(SCIP_PRESOL* presol)
{
   char name[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "unittest");

   if( strcmp(SCIPpresolGetName(presol), name) !=  0 )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

static
SCIP_RETCODE checkPresolGetDesc(SCIP_PRESOL* presol)
{
   char desc[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "presolver template");

   if( strcmp(SCIPpresolGetDesc(presol), desc) !=  0 )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

static
SCIP_RETCODE checkPresolGetPriority(SCIP_PRESOL* presol)
{
   if( SCIPpresolGetPriority(presol) !=  20010001 )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

static
SCIP_RETCODE checkPresolIsInitialized(SCIP_PRESOL* presol, SCIP_Bool val)
{
   if( SCIPpresolIsInitialized(presol) != val )
      return SCIP_ERROR;

   return SCIP_OKAY;
}


/** main function */
int main(
   int                        argc,
   char**                     argv
   )
{

   SCIP* scip;
   SCIP_RETCODE retcode;
   SCIP_PRESOL* presol;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-presol ===========\n");
   printf("=opt=  unittest-presol 0\n\n");

   /*********
    * Setup *
    *********/
   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /* include unittest presolver */
   SCIP_CALL( SCIPincludePresolUnittest(scip) );
   presol = SCIPfindPresol(scip, "unittest");

   /* create a problem and disable the presolver */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
   SCIP_CALL( initProb(scip) );

   /* set the msghdlr off */
   SCIPsetMessagehdlrQuiet(scip, TRUE);


   /*********
    * Tests *
    *********/

   /* tests before solving */
   CHECK_TEST( checkPresolGetName(presol) );
   CHECK_TEST( checkPresolGetDesc(presol) );
   CHECK_TEST( checkPresolGetPriority(presol) );
   CHECK_TEST( checkPresolIsInitialized(presol, FALSE) );

   /* solve */
   SCIP_CALL( SCIPsolve(scip) );
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   /* tests after solving */
   CHECK_TEST( checkPresolIsInitialized(presol, TRUE) );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return 0;
}
