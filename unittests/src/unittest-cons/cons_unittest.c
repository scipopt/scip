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

/**@file   cons_unittest.c
 * @brief  constraint handler for unittest constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_unittest.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "unittest"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP/**< propagation timing mask of the constraint handler*/

#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */

#define CONSHDLR_PRESOLTIMING  SCIP_PRESOLTIMING_ALWAYS


/* TODO: (optional) enable linear or nonlinear constraint upgrading */
#if 0
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#define LINCONSUPGD_PRIORITY          0 /**< priority of the constraint handler for upgrading of linear constraints */
#define NONLINCONSUPGD_PRIORITY       0 /**< priority of the constraint handler for upgrading of nonlinear constraints */
#endif


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for unittest constraints */
#if 0
struct SCIP_ConsData
{
};
#endif

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int nenfolp;   /* store the number of nenfolp calls */
   int ncheck;    /* store the number of check calls */
   int nsepalp;   /* store the number of sepalp calls */
   int nenfopslp; /* store the number of enfopslp calls */
   int nprop;     /* store the number of prop calls */
   int nresprop;  /* store the number of resprop calls */
   int npresol;   /* store the number of presol calls */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a unittest constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdUnittest)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to unittest constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to unittest constraint\n", SCIPconsGetName(cons));

      /* create the bin Unittest constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsUnittest(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}
#endif


/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyUnittest NULL
#endif

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


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitUnittest NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitUnittest NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreUnittest NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreUnittest NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolUnittest NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolUnittest NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteUnittest NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransUnittest NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpUnittest NULL
#endif


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



/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolUnittest NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpUnittest)
{

   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** vars;
   SCIP_ROW *row;
   SCIP_Bool infeasible;
   char s[SCIP_MAXSTRLEN];

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
   SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );
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


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveUnittest NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveUnittest NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableUnittest NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableUnittest NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsUnittest NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintUnittest NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyUnittest NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseUnittest NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest power constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsUnittest NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsUnittest)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of unittest power constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsUnittest NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for unittest constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create unittest constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   conshdlrdata->nenfolp = 0;
   conshdlrdata->ncheck = 0;
   conshdlrdata->nsepalp = 0;
   conshdlrdata->nenfopslp = 0;
   conshdlrdata->nprop = 0;
   conshdlrdata->nresprop = 0;
   conshdlrdata->npresol = 0;

   conshdlr = NULL;

   /* include constraint handler */
#if 0
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyUnittest,
         consFreeUnittest, consInitUnittest, consExitUnittest,
         consInitpreUnittest, consExitpreUnittest, consInitsolUnittest, consExitsolUnittest,
         consDeleteUnittest, consTransUnittest, consInitlpUnittest,
         consSepalpUnittest, consSepasolUnittest, consEnfolpUnittest, consEnfopsUnittest, consCheckUnittest,
         consPropUnittest, consPresolUnittest, consRespropUnittest, consLockUnittest,
         consActiveUnittest, consDeactiveUnittest,
         consEnableUnittest, consDisableUnittest, consDelvarsUnittest,
         consPrintUnittest, consCopyUnittest, consParseUnittest,
         consGetVarsUnittest, consGetNVarsUnittest, NULL, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpUnittest, consEnfopsUnittest, consCheckUnittest, consLockUnittest,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveUnittest) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyUnittest, consCopyUnittest) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveUnittest) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteUnittest) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsUnittest) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableUnittest) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableUnittest) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitUnittest) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreUnittest) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolUnittest) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeUnittest) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsUnittest) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsUnittest) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitUnittest) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreUnittest) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolUnittest) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpUnittest) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseUnittest) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolUnittest, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintUnittest) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropUnittest, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropUnittest) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpUnittest, consSepasolUnittest, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransUnittest) );
#endif

#ifdef LINCONSUPGD_PRIORITY
   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdUnittest, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
#endif

   /* add unittest constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a unittest constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsUnittest(
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsUnittest() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;


   /* find the unittest constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("unittest constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a unittest constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicUnittest(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            coefs,              /**< array with coefficients of constraint entries */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsUnittest(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/*
 * Interface methods of constraint handler
 */


/** gets nenfolp from the conshdlrdata */
int SCIPgetNenfolpUnittest(
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
int SCIPgetNcheckUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->ncheck;
}

/* gets nsepalp from the conshdlrdata */
int SCIPgetNsepalpUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nsepalp;
}

/* gets nenfopslp from the conshdlrdata */
int SCIPgetNenfopslpUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nenfopslp;
}

/* gets nprop from the conshdlrdata */
int SCIPgetNpropUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nprop;
}

/* gets nresprop from the conshdlrdata */
int SCIPgetNrespropUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->nresprop;
}

/* gets npresol from the conshdlrdata */
int SCIPgetNpresolUnittest(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   return conshdlrdata->npresol;
}
