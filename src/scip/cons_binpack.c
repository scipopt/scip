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

/**@file   cons_binpack.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for binpack constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_binpack.h"
#include "scip/cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "binpack"
#define CONSHDLR_DESC          "bin packing constraints of the form  a^T x <= b*y, x, y binary"
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

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyAnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   *valid = TRUE;

   return SCIP_OKAY;
}
#else
#define conshdlrCopyBinpack NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeBinpack NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitBinpack NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitBinpack NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreBinpack NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreBinpack NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolBinpack NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolBinpack NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteBinpack NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 0
static
SCIP_DECL_CONSTRANS(consTransBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransBinpack NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpBinpack NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpBinpack NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolBinpack NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropBinpack NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolBinpack NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropBinpack NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveBinpack NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveBinpack NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableBinpack NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableBinpack NULL
#endif

/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintBinpack)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintBinpack NULL
#endif

/** constraint copying method of constraint handler */
#define consCopyBinpack NULL

/** constraint parsing method of constraint handler */
#define consParseBinpack NULL


/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
static
SCIP_DECL_LINCONSUPGD(linconsUpgdBinpack)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a binpack constraint
    * - all variables must be binary
    * - all coefficients must be integral
    * - exactly one side is infinite
    * - there exists a coefficient a_k, s.t.
    *   - either left hand side is infinite, and right hand side is |a_k| + negcoeffsum,
    *     or right hand side is infinite, and left hand side is -|a_k| - negcoeffsum
    */
   upgrade = FALSE;
   if( nposbin + nnegbin == nvars
      && ncoeffspone + ncoeffsnone + ncoeffspint + ncoeffsnint == nvars
      && (SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs)) 
      && !(SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs)) )
   {
      SCIP_Bool found;
      int i;

      found = FALSE;
      if( SCIPisInfinity(scip, -lhs) )
      {
         assert(!SCIPisInfinity(scip, rhs));

         for( i = 0; i < nvars && !found; ++i )
            found = SCIPisEQ(scip, rhs, REALABS(vals[i]) + negcoeffsum);
      }
      else
      {
         assert(SCIPisInfinity(scip, rhs));

         for( i = 0; i < nvars && !found; ++i )
            found = SCIPisEQ(scip, lhs, -REALABS(vals[i]) - negcoeffsum);
      }

      upgrade = found;
   }

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to binpack constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Binpack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsBinpack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, rhs,
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

/** creates the handler for binpack constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBinpack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create binpack constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyBinpack,
         consFreeBinpack, consInitBinpack, consExitBinpack, 
         consInitpreBinpack, consExitpreBinpack, consInitsolBinpack, consExitsolBinpack,
         consDeleteBinpack, consTransBinpack, consInitlpBinpack,
         consSepalpBinpack, consSepasolBinpack, consEnfolpBinpack, consEnfopsBinpack, consCheckBinpack, 
         consPropBinpack, consPresolBinpack, consRespropBinpack, consLockBinpack,
         consActiveBinpack, consDeactiveBinpack, 
         consEnableBinpack, consDisableBinpack,
         consPrintBinpack, consCopyBinpack, consParseBinpack,
         conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdBinpack, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
#endif

   return SCIP_OKAY;
}

/** creates and captures a binpack constraint */
SCIP_RETCODE SCIPcreateConsBinpack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   len,                /**< number of nonzeros in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
   SCIP_Real*            vals,               /**< array with coefficients of constraint entries */
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
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of binpack constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the binpack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("binpack constraint handler not found\n");
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

