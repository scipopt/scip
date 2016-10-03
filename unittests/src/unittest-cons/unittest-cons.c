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

/* all methods in pub_cons.h
DONE:
SCIPconshdlrGetName
SCIPconshdlrGetDesc
SCIPconshdlrGetData
SCIPconshdlrGetSepaPriority
SCIPconshdlrGetEnfoPriority
SCIPconshdlrGetCheckPriority
SCIPconshdlrGetSepaFreq
SCIPconshdlrGetPropFreq
SCIPconshdlrGetEagerFreq
SCIPconshdlrNeedsCons
SCIPconshdlrDoesPresolve
SCIPconshdlrIsSeparationDelayed
SCIPconshdlrIsPropagationDelayed
SCIPconshdlrGetNEnfoLPCalls
SCIPconshdlrIsInitialized
SCIPconshdlrGetNCheckCalls
SCIPconshdlrGetNConss
SCIPconshdlrGetNEnfoConss
SCIPconshdlrGetNCheckConss
SCIPconshdlrGetNActiveConss
SCIPconshdlrGetNEnabledConss
SCIPconshdlrGetSetupTime
SCIPconshdlrGetPresolTime
SCIPconshdlrGetSepaTime
SCIPconshdlrGetEnfoLPTime
SCIPconshdlrGetEnfoPSTime
SCIPconshdlrGetPropTime
SCIPconshdlrGetStrongBranchPropTime
SCIPconshdlrGetCheckTime
SCIPconshdlrGetRespropTime
SCIPconshdlrGetNSepaCalls
SCIPconshdlrGetNEnfoPSCalls
SCIPconshdlrGetNPropCalls
SCIPconshdlrGetNRespropCalls
SCIPconshdlrGetNPresolCalls
SCIPconshdlrGetConss
SCIPconshdlrGetEnfoConss
SCIPconshdlrGetCheckConss


@TODO:
SCIPconshdlrGetNCutoffs
SCIPconshdlrGetNCutsFound
SCIPconshdlrGetNCutsApplied
SCIPconshdlrGetNConssFound
SCIPconshdlrGetNDomredsFound
SCIPconshdlrGetNChildren
SCIPconshdlrGetMaxNActiveConss
SCIPconshdlrGetStartNActiveConss
SCIPconshdlrGetNFixedVars
SCIPconshdlrGetNAggrVars
SCIPconshdlrGetNChgVarTypes
SCIPconshdlrGetNChgBds
SCIPconshdlrGetNAddHoles
SCIPconshdlrGetNDelConss
SCIPconshdlrGetNAddConss
SCIPconshdlrGetNUpgdConss
SCIPconshdlrGetNChgCoefs
SCIPconshdlrGetNChgSides
SCIPconshdlrWasLPSeparationDelayed
SCIPconshdlrWasSolSeparationDelayed
SCIPconshdlrWasPropagationDelayed
SCIPconshdlrIsClonable
SCIPconshdlrSetPropTiming
SCIPconshdlrGetPresolTiming
SCIPconshdlrSetPresolTiming
*/

static
SCIP_RETCODE consCheckName(SCIP_CONSHDLR* conshdlr )
{
   char name[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "unittest");

   CHECK_GET( strcmp(SCIPconshdlrGetName(conshdlr), name), 0 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckDesc(SCIP_CONSHDLR* conshdlr )
{
   char desc[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "constraint handler template");

   CHECK_GET( strcmp(SCIPconshdlrGetDesc(conshdlr), desc), 0 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckSepaPriority(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetSepaPriority(conshdlr), 0 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckEnfoPriority(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetEnfoPriority(conshdlr), 0 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckCheckPriority(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetCheckPriority(conshdlr), 0 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckSepaFreq(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetSepaFreq(conshdlr), -1 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckEagerFreq(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetEagerFreq(conshdlr), 100 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckPropFreq(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetPropFreq(conshdlr), -1 );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckNeedsCons(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrNeedsCons(conshdlr), TRUE );

   return SCIP_OKAY;
}


static
SCIP_RETCODE consCheckDoesPresolve(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrDoesPresolve(conshdlr), TRUE );

   return SCIP_OKAY;
}


static
SCIP_RETCODE consCheckIsSparationDelayed(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrIsSeparationDelayed(conshdlr), FALSE );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckIsPropagationDelayed(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrIsPropagationDelayed(conshdlr), FALSE );

   return SCIP_OKAY;
}

/*
static
SCIP_RETCODE consCheckWasLPSparationDelayed(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrWasLPSeparationDelayed(conshdlr),  );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckWasSolSparationDelayed(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrWasSolSeparationDelayed(conshdlr),  );
   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckWasPropagationDelayed(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrWasPropagationDelayed(conshdlr),  );
   return SCIP_OKAY;
}

*/

/*
static
SCIP_RETCODE consCheckIsCloneable(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrIsCloneable(conshdlr),  );
   return SCIP_OKAY;
}
*/

static
SCIP_RETCODE consCheckIsInitialized(SCIP_CONSHDLR* conshdlr, SCIP_Bool initialized)
{
   CHECK_GET( SCIPconshdlrIsInitialized(conshdlr), initialized );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetPropTiming(SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetPropTiming(conshdlr), SCIP_PROPTIMING_BEFORELP );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNConss(SCIP_CONSHDLR* conshdlr, int val)
{
   CHECK_GET( SCIPconshdlrGetNConss(conshdlr), val);

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNEnfoConss(SCIP_CONSHDLR* conshdlr, int val)
{
   CHECK_GET( SCIPconshdlrGetNEnfoConss(conshdlr), val);

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNCheckConss(SCIP_CONSHDLR* conshdlr, int val)
{
   CHECK_GET( SCIPconshdlrGetNCheckConss(conshdlr), val);

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNActiveConss(SCIP_CONSHDLR* conshdlr, int val)
{
   CHECK_GET( SCIPconshdlrGetNActiveConss(conshdlr), val);

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNEnabledConss(SCIP_CONSHDLR* conshdlr, int val)
{
   CHECK_GET( SCIPconshdlrGetNEnabledConss(conshdlr), val);

   return SCIP_OKAY;
}

/* how to test the time methods? */
static
SCIP_RETCODE consCheckGetSetupTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetSetupTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetPresolTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetPresolTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetSepaTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetSepaTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetEnfoLPTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetEnfoLPTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetEnfoPSTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetEnfoPSTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetPropTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetPropTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetStrongBranchPropTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetStrongBranchPropTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetCheckTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetCheckTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetRespropTime(SCIP_CONSHDLR* conshdlr)
{
   return (SCIPconshdlrGetRespropTime(conshdlr) >= 0) ? SCIP_OKAY : SCIP_ERROR;
}

static
SCIP_RETCODE consCheckGetNSepaCalls(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetNSepaCalls(conshdlr), SCIPgetNsepalpUnittest(scip) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetEnfoPSCalls(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetNEnfoPSCalls(conshdlr), SCIPgetNenfopslpUnittest(scip) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNPropCalls(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetNPropCalls(conshdlr), SCIPgetNpropUnittest(scip) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNRespropCalls(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetNRespropCalls(conshdlr), SCIPgetNrespropUnittest(scip) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckGetNPresolCalls(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetNPresolCalls(conshdlr), SCIPgetNpresolUnittest(scip) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE consCheckNEnfoLPCalls(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{

   CHECK_GET( SCIPconshdlrGetNEnfoLPCalls(conshdlr), SCIPgetNenfolpUnittest(scip) );

   return SCIP_OKAY;
}

#if 0
SCIP_RETCODE consCheckNCheckCalls(SCIP* scip, SCIP_CONSHDLR* conshdlr)
{
   CHECK_GET( SCIPconshdlrGetNCheckCalls(conshdlr), SCIPgetNcheckUnittest(scip) );

   return SCIP_OKAY;
}
#endif

/** main function */
int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP* scip;
   SCIP_RETCODE retcode;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS* cons;

   /* output stuff for automatic unittest evaluation */
   printf("@01 unittest-cons ===========\n");
   printf("=opt=  unittest-cons 0\n\n");

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
   conshdlr = SCIPfindConshdlr(scip, "unittest");


   /*********
    * Tests *
    *********/

   /* tests before solving */
   CHECK_TEST( consCheckName(conshdlr) );
   CHECK_TEST( consCheckDesc(conshdlr) );
   CHECK_TEST( consCheckSepaPriority(conshdlr) );
   CHECK_TEST( consCheckEnfoPriority(conshdlr) );
   CHECK_TEST( consCheckCheckPriority(conshdlr) );
   CHECK_TEST( consCheckSepaFreq(conshdlr) );
   CHECK_TEST( consCheckEagerFreq(conshdlr) );
   CHECK_TEST( consCheckNeedsCons(conshdlr) );
   CHECK_TEST( consCheckIsSparationDelayed(conshdlr) );
   CHECK_TEST( consCheckIsPropagationDelayed(conshdlr) );
   CHECK_TEST( consCheckGetPropTiming(conshdlr) );
   CHECK_TEST( consCheckIsInitialized(conshdlr, FALSE) );
   CHECK_TEST( consCheckDoesPresolve(conshdlr) );
   CHECK_TEST( consCheckPropFreq(conshdlr) );
   CHECK_TEST( consCheckGetNConss(conshdlr, 0) );
   CHECK_TEST( consCheckGetNEnfoConss(conshdlr, 0) );
   CHECK_TEST( consCheckGetNCheckConss(conshdlr, 0) );
   CHECK_TEST( consCheckGetNActiveConss(conshdlr, 0) );
   CHECK_TEST( consCheckGetNEnabledConss(conshdlr, 0) );
   CHECK_TEST( consCheckGetSetupTime(conshdlr) );
   CHECK_TEST( consCheckGetPresolTime(conshdlr) );
   CHECK_TEST( consCheckGetSepaTime(conshdlr) );
   CHECK_TEST( consCheckGetEnfoLPTime(conshdlr) );
   CHECK_TEST( consCheckGetEnfoPSTime(conshdlr) );
   CHECK_TEST( consCheckGetPropTime(conshdlr) );
   CHECK_TEST( consCheckGetEnfoLPTime(conshdlr) );
   CHECK_TEST( consCheckGetStrongBranchPropTime(conshdlr) );
   CHECK_TEST( consCheckGetCheckTime(conshdlr) );
   CHECK_TEST( consCheckGetRespropTime(conshdlr) );

   /* solve */
   SCIP_CALL( SCIPsolve(scip) );


   /* tests after solving */
   CHECK_TEST( consCheckIsInitialized(conshdlr, TRUE) );
   CHECK_TEST( consCheckNEnfoLPCalls(scip, conshdlr)  );
   CHECK_TEST( consCheckGetNActiveConss(conshdlr, 1)  );
   CHECK_TEST( consCheckGetNEnabledConss(conshdlr, 1)  );
   CHECK_TEST( consCheckGetSetupTime(conshdlr) );
   CHECK_TEST( consCheckGetPresolTime(conshdlr) );
   CHECK_TEST( consCheckGetSepaTime(conshdlr) );
   CHECK_TEST( consCheckGetEnfoLPTime(conshdlr) );
   CHECK_TEST( consCheckGetEnfoPSTime(conshdlr) );
   CHECK_TEST( consCheckGetPropTime(conshdlr) );
   CHECK_TEST( consCheckGetEnfoLPTime(conshdlr) );
   CHECK_TEST( consCheckGetStrongBranchPropTime(conshdlr) );
   CHECK_TEST( consCheckGetCheckTime(conshdlr) );
   CHECK_TEST( consCheckGetRespropTime(conshdlr) );
   CHECK_TEST( consCheckGetNSepaCalls(scip, conshdlr) );
   CHECK_TEST( consCheckGetEnfoPSCalls(scip, conshdlr) );
   CHECK_TEST( consCheckGetNPropCalls(scip, conshdlr) );
   CHECK_TEST( consCheckGetNRespropCalls(scip, conshdlr) );
   CHECK_TEST( consCheckGetNPresolCalls(scip, conshdlr) );
   CHECK_TEST( consCheckGetNEnfoConss(conshdlr, 1) );
   CHECK_TEST( consCheckGetNCheckConss(conshdlr, 1) );
   CHECK_TEST( consCheckGetNConss(conshdlr, 1) );

   /* We only count a call of the feasibility check method of a constraint handler if we check all constraints of a handler.
    * We want to compare this against SCIPconshdlrGetNCheckCalls(), but SCIP might call the check method of the constraint
    * handler to check a single constraint. In this case the counter for the number of check calls does not increase (for SCIP).
    * So the total number of calls of the check method should be at least SCIPconshdlrGetNCheckCalls(). */
   /*CHECK_TEST( consCheckNCheckCalls(scip, conshdlr) ); */


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
