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
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, 0.0, 2.0, FALSE, FALSE, TRUE) );
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

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

TestSuite(cons, .init = setup, .fini = teardown);
TestSuite(cons_solve, .init = setup_solve, .fini = teardown);


/* TESTS */

/* We only count a call of the feasibility check method of a constraint handler if we check all constraints of a handler.
 * We want to compare this against SCIPconshdlrGetNCheckCalls(), but SCIP might call the check method of the constraint
 * handler to check a single constraint. In this case the counter for the number of check calls does not increase (for SCIP).
 * So the total number of calls of the check method (getNcheckUnittests) should be at least SCIPconshdlrGetNCheckCalls().
 */
Test(cons, NCheckCalls)
{
   cr_assert_geq(getNcheckUnittest(scip), SCIPconshdlrGetNCheckCalls(conshdlr));
}

Test(cons, GetEnfoPriority)
{
   cr_assert_eq(SCIPconshdlrGetEnfoPriority(conshdlr), 0);
}

Test(cons, GetName)
{
   char name[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "unittest");
   cr_assert_str_eq(name, SCIPconshdlrGetName(conshdlr));
}

Test(cons, GetDesc)
{
   char desc[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(desc, SCIP_MAXSTRLEN, "constraint handler template");
   cr_assert_str_eq(desc, SCIPconshdlrGetDesc(conshdlr));
}

Test(cons, GetSepaPriority)
{
   cr_assert_eq(SCIPconshdlrGetSepaPriority(conshdlr), 0);
}

Test(cons, GetCheckPriority)
{
   cr_assert_eq(SCIPconshdlrGetCheckPriority(conshdlr), 0);
}

Test(cons, GetSepaFreq)
{
   cr_assert_eq(SCIPconshdlrGetSepaFreq(conshdlr), -1);
}

Test(cons, GetEagerFreq)
{
   cr_assert_eq(SCIPconshdlrGetEagerFreq(conshdlr), 100);
}

Test(cons, GetPropFreq)
{
   cr_assert_eq(SCIPconshdlrGetPropFreq(conshdlr), -1);
}

Test(cons, NeedsCons)
{
   cr_assert(SCIPconshdlrNeedsCons(conshdlr));
}

Test(cons, DoesPresolve)
{
   cr_assert(SCIPconshdlrDoesPresolve(conshdlr));
}

Test(cons, IsSeparationDelayed)
{
   cr_assert_not(SCIPconshdlrIsSeparationDelayed(conshdlr));
}

Test(cons, IsPropagationDelayed)
{
   cr_assert_not(SCIPconshdlrIsPropagationDelayed(conshdlr));
}

Test(cons, IsInitialized)
{
   cr_assert_not(SCIPconshdlrIsInitialized(conshdlr));
}

Test(cons, GetPropTiming)
{
   cr_assert_eq(SCIPconshdlrGetPropTiming(conshdlr), SCIP_PROPTIMING_BEFORELP);
}

Test(cons, GetNConss)
{
   cr_assert_eq(SCIPconshdlrGetNConss(conshdlr), 0);
}

Test(cons, GetNEnfoConss)
{
   cr_assert_eq(SCIPconshdlrGetNEnfoConss(conshdlr), 0);
}

Test(cons, GetNCheckConss)
{
   cr_assert_eq(SCIPconshdlrGetNCheckConss(conshdlr), 0);
}

Test(cons, GetNActiveConss)
{
   cr_assert_eq(SCIPconshdlrGetNActiveConss(conshdlr), 0);
}

Test(cons, GetNEnabledConss)
{
   cr_assert_eq(SCIPconshdlrGetNEnabledConss(conshdlr), 0);
}

Test(cons_solve, GetNEnabledConss)
{
   cr_assert_eq(SCIPconshdlrGetNEnabledConss(conshdlr), 1);
}

/* how to test the time methods? */
Test(cons_solve, GetSetupTime)
{
   cr_assert_geq(SCIPconshdlrGetSetupTime(conshdlr), 0.0);
}

Test(cons_solve, GetPresolTime)
{
   cr_assert_geq(SCIPconshdlrGetPresolTime(conshdlr), 0.0);
}

Test(cons_solve, GetSepaTime)
{
   cr_assert_geq(SCIPconshdlrGetSepaTime(conshdlr), 0.0);
}

Test(cons_solve, GetEnfoLPTime)
{
   cr_assert_geq(SCIPconshdlrGetEnfoLPTime(conshdlr), 0.0);
}

Test(cons_solve, GetEnfoPSTime)
{
   cr_assert_geq(SCIPconshdlrGetEnfoPSTime(conshdlr), 0.0);
}

Test(cons_solve, GetPropTime)
{
   cr_assert_geq(SCIPconshdlrGetPropTime(conshdlr), 0.0);
}

Test(cons, GetStrongBranchPropTime)
{
   cr_assert_geq(SCIPconshdlrGetStrongBranchPropTime(conshdlr), 0.0);
}

Test(cons, GetCheckTime)
{
   cr_assert_geq(SCIPconshdlrGetCheckTime(conshdlr), 0.0);
}

Test(cons, GetRespropTime)
{
   cr_assert_geq(SCIPconshdlrGetRespropTime(conshdlr), 0.0);
}

Test(cons, GetNSepaCalls)
{
   cr_assert_eq(SCIPconshdlrGetNSepaCalls(conshdlr), getNsepalpUnittest(scip));
}

Test(cons, GetEnfoPSCalls)
{
   cr_assert_eq(SCIPconshdlrGetNEnfoPSCalls(conshdlr), getNenfopslpUnittest(scip));
}

Test(cons, GetNPropCalls)
{
   cr_assert_eq(SCIPconshdlrGetNPropCalls(conshdlr), getNpropUnittest(scip));
}

Test(cons, GetNRespropCalls)
{
   cr_assert_eq(SCIPconshdlrGetNRespropCalls(conshdlr), getNrespropUnittest(scip));
}

Test(cons_solve, GetNPresolCalls)
{
   cr_assert_eq(SCIPconshdlrGetNPresolCalls(conshdlr), getNpresolUnittest(scip));
}

Test(cons_solve, NEnfoLPCalls)
{
   cr_assert_eq(SCIPconshdlrGetNEnfoLPCalls(conshdlr), getNenfolpUnittest(scip));
}

Test(cons_solve, IsInitialized)
{
   cr_assert(SCIPconshdlrIsInitialized(conshdlr));
}

Test(cons_solve, GetNConss)
{
   cr_assert_eq(SCIPconshdlrGetNConss(conshdlr), 1);
}

Test(cons_solve, GetNEnfoConss)
{
   cr_assert_eq(SCIPconshdlrGetNEnfoConss(conshdlr), 1);
}

Test(cons_solve, GetNCheckConss)
{
   cr_assert_eq(SCIPconshdlrGetNCheckConss(conshdlr), 1);
}

Test(cons_solve, GetNActiveConss)
{
   cr_assert_eq(SCIPconshdlrGetNActiveConss(conshdlr), 1);
}
