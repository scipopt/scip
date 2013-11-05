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

/**@file   unittest-cons.c
 * @brief  unit test for checking setters on scip.c
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include "cons_unittest.h"
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


/* local methods */

/** create bounded problem */
static
SCIP_RETCODE initProb(
   SCIP*                 scip                /**< SCIP instance */
   )
{
   SCIP_CONS* cons;
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
static
SCIP_RETCODE consCheckName(SCIP* scip, SCIP_CONSHDLR* hdler )
{
   const char* name = "unittest";
   CHECK_GET( strcmp(SCIPconshdlrGetName(hdler), name), 0 );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckDesc(SCIP* scip, SCIP_CONSHDLR* hdler )
{
   const char* desc = "constraint handler template";
   CHECK_GET( strcmp(SCIPconshdlrGetDesc(hdler), desc), 0 );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckSepaPriority(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrGetSepaPriority(hdler), 0 );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckEnfoPriority(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrGetEnfoPriority(hdler), 0 );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckCheckPriority(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrGetCheckPriority(hdler), 0 );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckSepaFreq(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrGetSepaFreq(hdler), -1 );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckEagerFreq(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrGetEagerFreq(hdler), 100 );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckNeedsCons(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrNeedsCons(hdler), TRUE );
   return SCIP_OKAY;
}

/*
static
SCIP_RETCODE consCheckDoesPresolve(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrDoesPresolve(hdler),  );
   return SCIP_OKAY;
}
*/

static
SCIP_RETCODE consCheckIsSparationDelayed(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrIsSeparationDelayed(hdler), FALSE );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckIsPropagationDelayed(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrIsPropagationDelayed(hdler), FALSE );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckIsPresolvingDelayed(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrIsPresolvingDelayed(hdler), FALSE );
   return SCIP_OKAY;
}

/*
static
SCIP_RETCODE consCheckWasLPSparationDelayed(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrWasLPSeparationDelayed(hdler),  );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckWasSolSparationDelayed(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrWasSolSeparationDelayed(hdler),  );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckWasPropagationDelayed(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrWasPropagationDelayed(hdler),  );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckWasPresolvingDelayed(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrWasPresolvingDelayed(hdler),  );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckIsInitialized(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrIsInitialized(hdler),  );
   return SCIP_OKAY;
}
*/

/*
static
SCIP_RETCODE consCheckIsCloneable(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrIsCloneable(hdler),  );
   return SCIP_OKAY;
}
*/

static
SCIP_RETCODE consCheckGetPropTimingmask(SCIP* scip, SCIP_CONSHDLR* hdler)
{
   CHECK_GET( SCIPconshdlrGetPropTimingmask(hdler), SCIP_PROPTIMING_BEFORELP );
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
   SCIP* scip;
   SCIP_CONSHDLR* hdler;
   SCIP_CONS* cons;

   /*********
    * Setup *
    *********/
   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include unittest constraint handler */
   SCIP_CALL( SCIPincludeConshdlrUnittest(scip) );

   /* create a problem and disable the presolver */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
   SCIP_CALL( initProb(scip) );

   /* set the msghdlr off */
   SCIPsetMessagehdlrQuiet(scip, TRUE);

   /* create a constraint of the unittesthandler */
   /* the constraint handler just add the constraint: x+y <=2 */
   SCIP_CALL( SCIPcreateConsUnittest(scip, &cons, "UC", 2, NULL, NULL, 0, 2, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
         FALSE, FALSE, FALSE, FALSE));

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* get the constraint handler */
   hdler = SCIPfindConshdlr(scip,"unittest");


   /*********
    * Tests *
    *********/

   /* tests before solving */
   consCheckName(scip, hdler);
   consCheckDesc(scip, hdler);
   consCheckSepaPriority(scip, hdler);
   consCheckEnfoPriority(scip, hdler);
   consCheckCheckPriority(scip, hdler);
   consCheckSepaFreq(scip, hdler);
   consCheckEagerFreq(scip, hdler);
   consCheckNeedsCons(scip, hdler);
   consCheckIsSparationDelayed(scip, hdler);
   consCheckIsPropagationDelayed(scip, hdler);
   consCheckIsPresolvingDelayed(scip, hdler);
   consCheckGetPropTimingmask(scip, hdler);

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
