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

/**@file   cons.c
 * @brief  unit test for checking setters on scip.c
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include "scip/scipdefplugins.h"
#include <string.h>

/* UNIT TEST CONSHDLR */
/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "unittest"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */

#define CONSHDLR_PRESOLTIMING  SCIP_PRESOLTIMING_ALWAYS

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   nenfolp;            /**< store the number of nenfolp calls */
   int                   ncheck;             /**< store the number of check calls */
   int                   nsepalp;            /**< store the number of sepalp calls */
   int                   nenfopslp;          /**< store the number of enfopslp calls */
   int                   nprop;              /**< store the number of prop calls */
   int                   nresprop;           /**< store the number of resprop calls */
   int                   npresol;            /**< store the number of presol calls */
};


/*
 * Callback methods of constraint handler
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeUnittest)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpUnittest)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->nsepalp++;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpUnittest)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** vars;
   SCIP_ROW *row;
   SCIP_Bool infeasible;
   char s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( result != NULL );

   /* count this function call in the hdlr data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   assert(conshdlrdata != NULL);

   conshdlrdata->nenfolp++;

   /* now add a cutting plane: x+y <= 2 */
   vars = SCIPgetVars(scip);
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "MyCut");
   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, s, 0.0, 2.0, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[0], 1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, vars[0], 1.0) );
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );
   SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
   SCIP_CALL( SCIPreleaseRow(scip, &row));

   *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsUnittest)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->nenfopslp++;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckUnittest)
{
   SCIP_Real val;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** vars;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( result != NULL );

   vars = SCIPgetVars(scip);

   assert(vars != NULL);
   assert(vars[0] != NULL);
   assert(vars[1] != NULL);

   /* count this function call in the hdlr data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   assert(conshdlrdata != NULL);

   conshdlrdata->ncheck++;

   val = SCIPgetSolVal(scip, sol, vars[0]) + SCIPgetSolVal(scip, sol, vars[1]);

   if( val > 2 )
      *result = SCIP_INFEASIBLE;
   else
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropUnittest)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->nprop++;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolUnittest)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->npresol++;
   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropUnittest)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->nresprop++;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockUnittest)
{
   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for unittest constraints and includes it in SCIP */
static
SCIP_RETCODE includeConshdlrUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr = NULL;

   /* create unittest constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   conshdlrdata->nenfolp = 0;
   conshdlrdata->ncheck = 0;
   conshdlrdata->nsepalp = 0;
   conshdlrdata->nenfopslp = 0;
   conshdlrdata->nprop = 0;
   conshdlrdata->nresprop = 0;
   conshdlrdata->npresol = 0;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpUnittest, consEnfopsUnittest, consCheckUnittest, consLockUnittest,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeUnittest) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolUnittest, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropUnittest, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropUnittest) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpUnittest, NULL, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );

   return SCIP_OKAY;
}

/** creates and captures a unittest constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
static
SCIP_RETCODE createConsUnittest(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata = NULL;

   /* find the unittest constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("unittest constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}


/*
 * Interface methods of constraint handler
 */

/** gets nenfolp from the conshdlrdata */
static
int getNenfolpUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nenfolp;
}

/** gets nenfolp from the conshdlrdata */
static
int getNcheckUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->ncheck;
}

/** gets nsepalp from the conshdlrdata */
static
int getNsepalpUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nsepalp;
}

/** gets nenfopslp from the conshdlrdata */
static
int getNenfopslpUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nenfopslp;
}

/** gets nprop from the conshdlrdata */
static
int getNpropUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nprop;
}

/** gets nresprop from the conshdlrdata */
static
int getNrespropUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nresprop;
}

/** gets npresol from the conshdlrdata */
static
int getNpresolUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->npresol;
}

/** END CONSHDLR **/

#include "include/scip_test.h"

/*
 * HELPER METHODS
 */


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

/* GLOBAL VARIABLES */
static SCIP_CONSHDLR* conshdlr;
static SCIP* scip;

/* TEST SUITES */

/** setup of test run */
static
void setup(void)
{
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_CONS* cons;

   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include unittest constraint handler */
   SCIP_CALL( includeConshdlrUnittest(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0, 2, -1.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0, 2, -1.0, SCIP_VARTYPE_INTEGER) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );

   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );

   /* create a constraint of the unittesthandler: it just adds the constraint x + y <= 2 */
   SCIP_CALL( createConsUnittest(scip, &cons, "UC", 2, NULL, NULL, 0, 2, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
         FALSE, FALSE, FALSE, FALSE));

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* set the msghdlr off */
   SCIPsetMessagehdlrQuiet(scip, TRUE);

   /* get the constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "unittest");
}

/** setup solving test */
static
void setup_solve(void)
{
   setup();

   /* solve */
   SCIP_CALL( SCIPsolve(scip) );
}

/** deinitialization method */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}



/* TESTS */

/* We only count a call of the feasibility check method of a constraint handler if we check all constraints of a handler.
 * We want to compare this against SCIPconshdlrGetNCheckCalls(), but SCIP might call the check method of the constraint
 * handler to check a single constraint. In this case the counter for the number of check calls does not increase (for SCIP).
 * So the total number of calls of the check method (getNcheckUnittests) should be at least SCIPconshdlrGetNCheckCalls().
 */
void test_cons_NCheckCalls(void)
{
   TEST_ASSERT_GREATER_OR_EQUAL(getNcheckUnittest(scip), SCIPconshdlrGetNCheckCalls(conshdlr));
}

void test_cons_GetEnfoPriority(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetEnfoPriority(conshdlr), 0);
}

void test_cons_GetName(void)
{
   char name[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "unittest");
   TEST_ASSERT_EQUAL_STRING(name, SCIPconshdlrGetName(conshdlr));
}

void test_cons_GetDesc(void)
{
   char desc[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "constraint handler template");
   TEST_ASSERT_EQUAL_STRING(desc, SCIPconshdlrGetDesc(conshdlr));
}

void test_cons_GetSepaPriority(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetSepaPriority(conshdlr), 0);
}

void test_cons_GetCheckPriority(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetCheckPriority(conshdlr), 0);
}

void test_cons_GetSepaFreq(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetSepaFreq(conshdlr), -1);
}

void test_cons_GetEagerFreq(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetEagerFreq(conshdlr), 100);
}

void test_cons_GetPropFreq(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetPropFreq(conshdlr), -1);
}

void test_cons_NeedsCons(void)
{
   TEST_ASSERT(SCIPconshdlrNeedsCons(conshdlr));
}

void test_cons_DoesPresolve(void)
{
   TEST_ASSERT(SCIPconshdlrDoesPresolve(conshdlr));
}

void test_cons_IsSeparationDelayed(void)
{
   TEST_ASSERT_NOT(SCIPconshdlrIsSeparationDelayed(conshdlr));
}

void test_cons_IsPropagationDelayed(void)
{
   TEST_ASSERT_NOT(SCIPconshdlrIsPropagationDelayed(conshdlr));
}

void test_cons_IsInitialized(void)
{
   TEST_ASSERT_NOT(SCIPconshdlrIsInitialized(conshdlr));
}

void test_cons_GetPropTiming(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetPropTiming(conshdlr), SCIP_PROPTIMING_BEFORELP);
}

void test_cons_GetNConss(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNConss(conshdlr), 0);
}

void test_cons_GetNEnfoConss(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNEnfoConss(conshdlr), 0);
}

void test_cons_GetNCheckConss(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNCheckConss(conshdlr), 0);
}

void test_cons_GetNActiveConss(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNActiveConss(conshdlr), 0);
}

void test_cons_GetNEnabledConss(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNEnabledConss(conshdlr), 0);
}

/* Forward declaration for solve tests */
static void ensureSolved(void);

void test_cons_solve_GetNEnabledConss(void)
{
   ensureSolved();
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNEnabledConss(conshdlr), 1);
}

/* how to test the time methods? */
void test_cons_solve_GetSetupTime(void)
{
   ensureSolved();
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetSetupTime(conshdlr), 0.0);
}

void test_cons_solve_GetPresolTime(void)
{
   ensureSolved();
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetPresolTime(conshdlr), 0.0);
}

void test_cons_solve_GetSepaTime(void)
{
   ensureSolved();
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetSepaTime(conshdlr), 0.0);
}

void test_cons_solve_GetEnfoLPTime(void)
{
   ensureSolved();
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetEnfoLPTime(conshdlr), 0.0);
}

void test_cons_solve_GetEnfoPSTime(void)
{
   ensureSolved();
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetEnfoPSTime(conshdlr), 0.0);
}

void test_cons_solve_GetPropTime(void)
{
   ensureSolved();
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetPropTime(conshdlr), 0.0);
}

void test_cons_GetStrongBranchPropTime(void)
{
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetStrongBranchPropTime(conshdlr), 0.0);
}

void test_cons_GetCheckTime(void)
{
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetCheckTime(conshdlr), 0.0);
}

void test_cons_GetRespropTime(void)
{
   TEST_ASSERT_GREATER_OR_EQUAL(SCIPconshdlrGetRespropTime(conshdlr), 0.0);
}

void test_cons_GetNSepaCalls(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNSepaCalls(conshdlr), getNsepalpUnittest(scip));
}

void test_cons_GetEnfoPSCalls(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNEnfoPSCalls(conshdlr), getNenfopslpUnittest(scip));
}

void test_cons_GetNPropCalls(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNPropCalls(conshdlr), getNpropUnittest(scip));
}

void test_cons_GetNRespropCalls(void)
{
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNRespropCalls(conshdlr), getNrespropUnittest(scip));
}

void test_cons_solve_GetNPresolCalls(void)
{
   ensureSolved();
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNPresolCalls(conshdlr), getNpresolUnittest(scip));
}

void test_cons_solve_NEnfoLPCalls(void)
{
   ensureSolved();
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNEnfoLPCalls(conshdlr), getNenfolpUnittest(scip));
}

/* Helper: ensure SCIP is solved for solve tests */
static void ensureSolved(void)
{
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVED )
   {
      SCIP_CALL( SCIPsolve(scip) );
   }
}

void test_cons_solve_IsInitialized(void)
{
   ensureSolved();
   TEST_ASSERT(SCIPconshdlrIsInitialized(conshdlr));
}

void test_cons_solve_GetNConss(void)
{
   ensureSolved();
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNConss(conshdlr), 1);
}

void test_cons_solve_GetNEnfoConss(void)
{
   ensureSolved();
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNEnfoConss(conshdlr), 1);
}

void test_cons_solve_GetNCheckConss(void)
{
   ensureSolved();
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNCheckConss(conshdlr), 1);
}

void test_cons_solve_GetNActiveConss(void)
{
   ensureSolved();
   TEST_ASSERT_EQUAL(SCIPconshdlrGetNActiveConss(conshdlr), 1);
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_cons_NCheckCalls);
   RUN_TEST(test_cons_GetEnfoPriority);
   RUN_TEST(test_cons_GetName);
   RUN_TEST(test_cons_GetDesc);
   RUN_TEST(test_cons_GetSepaPriority);
   RUN_TEST(test_cons_GetCheckPriority);
   RUN_TEST(test_cons_GetSepaFreq);
   RUN_TEST(test_cons_GetEagerFreq);
   RUN_TEST(test_cons_GetPropFreq);
   RUN_TEST(test_cons_NeedsCons);
   RUN_TEST(test_cons_DoesPresolve);
   RUN_TEST(test_cons_IsSeparationDelayed);
   RUN_TEST(test_cons_IsPropagationDelayed);
   RUN_TEST(test_cons_IsInitialized);
   RUN_TEST(test_cons_GetPropTiming);
   RUN_TEST(test_cons_GetNConss);
   RUN_TEST(test_cons_GetNEnfoConss);
   RUN_TEST(test_cons_GetNCheckConss);
   RUN_TEST(test_cons_GetNActiveConss);
   RUN_TEST(test_cons_GetNEnabledConss);
   RUN_TEST(test_cons_solve_GetNEnabledConss);
   RUN_TEST(test_cons_solve_GetSetupTime);
   RUN_TEST(test_cons_solve_GetPresolTime);
   RUN_TEST(test_cons_solve_GetSepaTime);
   RUN_TEST(test_cons_solve_GetEnfoLPTime);
   RUN_TEST(test_cons_solve_GetEnfoPSTime);
   RUN_TEST(test_cons_solve_GetPropTime);
   RUN_TEST(test_cons_GetStrongBranchPropTime);
   RUN_TEST(test_cons_GetCheckTime);
   RUN_TEST(test_cons_GetRespropTime);
   RUN_TEST(test_cons_GetNSepaCalls);
   RUN_TEST(test_cons_GetEnfoPSCalls);
   RUN_TEST(test_cons_GetNPropCalls);
   RUN_TEST(test_cons_GetNRespropCalls);
   RUN_TEST(test_cons_solve_GetNPresolCalls);
   RUN_TEST(test_cons_solve_NEnfoLPCalls);
   RUN_TEST(test_cons_solve_IsInitialized);
   RUN_TEST(test_cons_solve_GetNConss);
   RUN_TEST(test_cons_solve_GetNEnfoConss);
   RUN_TEST(test_cons_solve_GetNCheckConss);
   RUN_TEST(test_cons_solve_GetNActiveConss);
   return UNITY_END();
}
