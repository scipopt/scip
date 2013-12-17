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

/**@file   unittest-presol.c
 * @brief  unit test for checking setters on scip.c
 * @author Benjamin Mueller
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



/* local methods */

/** create bounded problem */
static
SCIP_RETCODE initProb(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0, 2, -1.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0, 2, -1.0, SCIP_VARTYPE_INTEGER) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );

   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );

   return SCIP_OKAY;
}


/* Check methods */




/************************************/

/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{

   SCIP* scip;
   SCIP_RETCODE retcode;

   /*********
    * Setup *
    *********/
   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /* create a problem and disable the presolver */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
   SCIP_CALL( initProb(scip) );

   /* set the msghdlr off */
   SCIPsetMessagehdlrQuiet(scip, TRUE);




   /*********
    * Tests *
    *********/

   /* tests before solving */


   /* solve */
   SCIP_CALL( SCIPsolve(scip) );



   /* tests after solving */



   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   printf("All tests passed\n");
   return 0;
}
