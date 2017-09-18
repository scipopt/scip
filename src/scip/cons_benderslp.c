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

/**@file   cons_benderslp.c
 * @brief  constraint handler for benderslp decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/scip.h"

#include "scip/cons_benderslp.h"
#include "scip/cons_benders.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "benderslp"
#define CONSHDLR_DESC          "constraint handler for Benders' Decomposition to separate root node LP solutions"
#define CONSHDLR_ENFOPRIORITY     10000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY    10000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for benderslp constraints */
//struct SCIP_ConsData
//{
//};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   ncalls;             /**< the number of calls to the constraint handler. */
};


/*
 * Local methods
 */

#if 0
/* applies the generated cut to the master problem*/
static
SCIP_RETCODE applyCut(
   //SCIP_BENDERSLP_CUTTYPE  cuttype
   )
{
   return SCIP_OKAY;
}

/* apply Magnanti-Wong strengthening of the dual solutions from the optimal LP */
static
SCIP_RETCODE applyMagnantiWongDualStrengthening(
   )
{
   return SCIP_OKAY;
}
#endif


/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 1
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBenderslp)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeConshdlrBenderslp(scip) );

   return SCIP_OKAY;
}
#else
#define conshdlrCopyBenderslp NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBenderslp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* free memory for conshdlrdata*/
   if( conshdlrdata != NULL )
   {
      SCIPfreeMemory(scip, &conshdlrdata);
   }

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   /* A possible implementation from gcg_pricer.cpp */
   assert(scip == scip_);
   assert(reducedcostpricing != NULL);
   assert(farkaspricing != NULL);

   SCIP_CALL( solversInit() );

   SCIP_CALL( reducedcostpricing->resetCalls() );
   SCIP_CALL( farkaspricing->resetCalls() );

   return SCIP_OKAY;
}
#else
#define consInitBenderslp NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitBenderslp NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreBenderslp NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreBenderslp NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolBenderslp NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolBenderslp NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteBenderslp NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransBenderslp NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpBenderslp NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpBenderslp NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolBenderslp NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBenderslp)
{  /*lint --e{715}*/

   if( SCIPgetDepth(scip) >= 1 )
      (*result) = SCIP_FEASIBLE;
   else
      SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, NULL, conshdlr, result, LP) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBenderslp)
{  /*lint --e{715}*/

   if( SCIPgetDepth(scip) >= 1 )
      (*result) = SCIP_FEASIBLE;
   else
      SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, sol, conshdlr, result, RELAX) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBenderslp)
{  /*lint --e{715}*/

   if( SCIPgetDepth(scip) >= 1 )
      (*result) = SCIP_FEASIBLE;
   else
      SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, NULL, conshdlr, result, PSEUDO) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
/*  This function checks the feasibility of the Benderslp' decomposition master problem. In the case that the problem is
 *  feasible, then the auxiliary variables must be updated with the subproblem objective function values. The update
 *  occurs in the solve subproblems function. */
static
SCIP_DECL_CONSCHECK(consCheckBenderslp)
{  /*lint --e{715}*/
   (*result) = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropBenderslp NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolBenderslp NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropBenderslp NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBenderslp)
{  /*lint --e{715}*/
   //SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   //SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveBenderslp NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveBenderslp NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableBenderslp NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableBenderslp NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsBenderslp NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintBenderslp NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyBenderslp NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseBenderslp NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsBenderslp NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsBenderslp NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsBenderslp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsBenderslp NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for benderslp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBenderslp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   /* create benderslp constraint handler data */
   conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->ncalls = 0;

   conshdlr = NULL;

   /* include constraint handler */
#if 0
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyBenderslp,
         consFreeBenderslp, consInitBenderslp, consExitBenderslp,
         consInitpreBenderslp, consExitpreBenderslp, consInitsolBenderslp, consExitsolBenderslp,
         consDeleteBenderslp, consTransBenderslp, consInitlpBenderslp,
         consSepalpBenderslp, consSepasolBenderslp, consEnfolpBenderslp, consEnforelaxBenderslp, consEnfopsBenderslp, consCheckBenderslp,
         consPropBenderslp, consPresolBenderslp, consRespropBenderslp, consLockBenderslp,
         consActiveBenderslp, consDeactiveBenderslp,
         consEnableBenderslp, consDisableBenderslp, consDelvarsBenderslp,
         consPrintBenderslp, consCopyBenderslp, consParseBenderslp,
         consGetVarsBenderslp, consGetNVarsBenderslp, consGetDiveBdChgsBenderslp, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBenderslp, consEnfopsBenderslp, consCheckBenderslp, consLockBenderslp,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveBenderslp) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBenderslp, consCopyBenderslp) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveBenderslp) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteBenderslp) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsBenderslp) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableBenderslp) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableBenderslp) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitBenderslp) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreBenderslp) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolBenderslp) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBenderslp) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsBenderslp) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsBenderslp) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsBenderslp) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitBenderslp) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreBenderslp) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolBenderslp) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpBenderslp) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseBenderslp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolBenderslp, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintBenderslp) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropBenderslp, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropBenderslp) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpBenderslp, consSepasolBenderslp, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransBenderslp) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBenderslp) );

#endif


   return SCIP_OKAY;
}

/** creates and captures a benderslp constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBenderslp(
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsBenderslp() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of benderslp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the benderslp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("benderslp constraint handler not found\n");
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

/** creates and captures a benderslp constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicBenderslp(
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
   SCIP_CALL( SCIPcreateConsBenderslp(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
