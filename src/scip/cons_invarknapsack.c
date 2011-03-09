/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_invarknapsack.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for invarknapsack constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_invarknapsack.h"
#include "scip/cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "invarknapsack"
#define CONSHDLR_DESC          "invariant knapsack constraint of the form  1^T x <= b or 1^T x == b, x binary"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define LINCONSUPGD_PRIORITY    +000000 /**< priority of the constraint handler for upgrading of linear constraints */




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   *valid = FALSE;

   return SCIP_OKAY;
}
#else
#define conshdlrCopyInvarknapsack NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeInvarknapsack NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitInvarknapsack NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitInvarknapsack NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreInvarknapsack NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreInvarknapsack NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolInvarknapsack NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolInvarknapsack NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteInvarknapsack NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 0
static
SCIP_DECL_CONSTRANS(consTransInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransInvarknapsack NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpInvarknapsack NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpInvarknapsack NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolInvarknapsack NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropInvarknapsack NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolInvarknapsack NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropInvarknapsack NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveInvarknapsack NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveInvarknapsack NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableInvarknapsack NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableInvarknapsack NULL
#endif

/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintInvarknapsack NULL
#endif

/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyInvarknapsack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   *valid = FALSE;

   return SCIP_OKAY;
}
#else
#define consCopyInvarknapsack NULL
#endif


/** constraint parsing method of constraint handler */
#define consParseInvarknapsack NULL


/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
static
SCIP_DECL_LINCONSUPGD(linconsUpgdInvarknapsack)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a invariant knapsack constraint
    * - all coefficients must be +/- 1
    * - all variables must be binary
    * - either one of the sides is infinite, or both sides are equal
    */
   upgrade = (nposbin + nnegbin == nvars)
      && (ncoeffspone + ncoeffsnone == nvars)
      && (SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs) || SCIPisEQ(scip, lhs, rhs));

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to invarknapsack constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Invarknapsack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsInvarknapsack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}
#endif




/*
 * constraint specific interface methods
 */

/** creates the handler for invarknapsack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrInvarknapsack(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create invarknapsack constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         conshdlrCopyInvarknapsack,
         consFreeInvarknapsack, consInitInvarknapsack, consExitInvarknapsack, 
         consInitpreInvarknapsack, consExitpreInvarknapsack, consInitsolInvarknapsack, consExitsolInvarknapsack,
         consDeleteInvarknapsack, consTransInvarknapsack, consInitlpInvarknapsack,
         consSepalpInvarknapsack, consSepasolInvarknapsack, consEnfolpInvarknapsack, consEnfopsInvarknapsack, 
         consCheckInvarknapsack, 
         consPropInvarknapsack, consPresolInvarknapsack, consRespropInvarknapsack, consLockInvarknapsack,
         consActiveInvarknapsack, consDeactiveInvarknapsack, 
         consEnableInvarknapsack, consDisableInvarknapsack,
         consPrintInvarknapsack, consCopyInvarknapsack, consParseInvarknapsack,
         conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdInvarknapsack, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
#endif

   return SCIP_OKAY;
}

/** creates and captures a invarknapsack constraint */
SCIP_RETCODE SCIPcreateConsInvarknapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   len,                /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
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
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of invarknapsack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the invarknapsack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("invarknapsack constraint handler not found\n");
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

