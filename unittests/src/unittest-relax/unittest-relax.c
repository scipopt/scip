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

/**@file   unittest-relax.c
 * @brief  unit test for checking setters on scip.c
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include "scip/scipdefplugins.h"
#include "relax_unittest.h"
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

/** macro to check the value of a 'getter' and 'value'
  *
  */
#define CHECK_GET(getter, value)   \
   do                              \
   {                               \
      if( getter != value )        \
         return SCIP_ERROR;        \
   }                               \
   while(FALSE)




/* Check methods */
static
SCIP_RETCODE relaxCheckName(SCIP_RELAX* relax)
{
   char name[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "relax-unittest");

   CHECK_GET( strcmp(SCIPrelaxGetName(relax), name), 0 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE relaxCheckDesc(SCIP_RELAX* relax)
{
   char desc[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "relaxator template");

   CHECK_GET( strcmp(SCIPrelaxGetDesc(relax), desc), 0 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE relaxCheckPriority(SCIP_RELAX* relax)
{
   int priority;

   priority = 101;

   CHECK_GET( SCIPrelaxGetPriority(relax), priority );

   return SCIP_OKAY;
}

static
SCIP_RETCODE relaxCheckFreq(SCIP_RELAX* relax)
{
   int freq;

   freq = 2;

   CHECK_GET( SCIPrelaxGetFreq(relax), freq );

   return SCIP_OKAY;
}

/*@todo how to check this? */
static
SCIP_RETCODE relaxCheckSetupTime(SCIP_RELAX* relax)
{
   if( SCIPrelaxGetSetupTime(relax) < 0 )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

/*@todo how to check this? */
static
SCIP_RETCODE relaxCheckTime(SCIP_RELAX* relax)
{
   if( SCIPrelaxGetTime(relax) < 0 )
      return SCIP_ERROR;

   return SCIP_OKAY;
}

/*@todo how to check this? */
static
SCIP_RETCODE relaxCheckNCalls(SCIP* scip, SCIP_RELAX* relax)
{
   CHECK_GET( SCIPrelaxGetNCalls(relax), SCIPgetNcallsUnittest(scip) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE relaxCheckInitialized(SCIP_RELAX* relax, SCIP_Bool val)
{
   CHECK_GET( SCIPrelaxIsInitialized(relax), val );

   return SCIP_OKAY;
}

/*@todo how to check this? */
static
SCIP_RETCODE relaxCheckMarkedUnsolved(SCIP_RELAX* relax)
{
   return SCIP_OKAY;
}

/** main function */
int main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;
   SCIP* scip;
   SCIP_RELAX* relax;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-relax ===========\n");
   printf("=opt=  unittest-relax 0\n\n");

   /*********
    * Setup *
    *********/
   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include binpacking reader */
   SCIP_CALL( SCIPincludeRelaxUnittest(scip) );

   /* store the relaxertor */
   relax = SCIPfindRelax(scip, "relax-unittest");

   /* create a problem and disable the presolver */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* set the msghdlr off */
   SCIPsetMessagehdlrQuiet(scip, TRUE);

   /*********
    * Tests *
    *********/
   CHECK_TEST( relaxCheckName(relax) );
   CHECK_TEST( relaxCheckDesc(relax) );
   CHECK_TEST( relaxCheckFreq(relax) );
   CHECK_TEST( relaxCheckTime(relax) );
   CHECK_TEST( relaxCheckPriority(relax) );
   CHECK_TEST( relaxCheckSetupTime(relax) );
   CHECK_TEST( relaxCheckNCalls(scip, relax) );

   CHECK_TEST( relaxCheckInitialized(relax, FALSE) );

   SCIP_CALL( SCIPsolve(scip) );

   CHECK_TEST( relaxCheckInitialized(relax, TRUE) );
   CHECK_TEST( relaxCheckMarkedUnsolved(relax) );

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
