/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_knapsack.c
 * @brief  constraint handler for knapsack constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_knapsack.h"
#include "cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "knapsack"
#define CONSHDLR_DESC          "knapsack constraint of the form  a^T x <= b, x binary"
#define CONSHDLR_SEPAPRIORITY   +000000
#define CONSHDLR_ENFOPRIORITY   +000000
#define CONSHDLR_CHECKPRIORITY  +000000
#define CONSHDLR_SEPAFREQ            -1
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE

#define LINCONSUPGD_PRIORITY    +000000




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
DECL_CONSFREE(consFreeKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeKnapsack NULL
#endif


/** initialization method of constraint handler (called when problem solving starts) */
#if 0
static
DECL_CONSINIT(consInitKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitKnapsack NULL
#endif


/** deinitialization method of constraint handler (called when problem solving exits) */
#if 0
static
DECL_CONSEXIT(consExitKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitKnapsack NULL
#endif


/** solving start notification method of constraint handler (called when presolving was finished) */
#if 0
static
DECL_CONSSOLSTART(consSolstartKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSolstartKnapsack NULL
#endif


/** frees specific constraint data */
#if 0
static
DECL_CONSDELETE(consDeleteKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteKnapsack NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 0
static
DECL_CONSTRANS(consTransKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransKnapsack NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
DECL_CONSINITLP(consInitlpKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpKnapsack NULL
#endif


/** separation method of constraint handler */
#if 0
static
DECL_CONSSEPA(consSepaKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepaKnapsack NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
DECL_CONSPROP(consPropKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropKnapsack NULL
#endif


/** presolving method of constraint handler */
#if 0
static
DECL_CONSPRESOL(consPresolKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolKnapsack NULL
#endif


/** conflict variable resolving method of constraint handler */
#if 0
static
DECL_CONSRESCVAR(consRescvarKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRescvarKnapsack NULL
#endif


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
DECL_CONSACTIVE(consActiveKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveKnapsack NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
DECL_CONSDEACTIVE(consDeactiveKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveKnapsack NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
DECL_CONSENABLE(consEnableKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableKnapsack NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
DECL_CONSDISABLE(consDisableKnapsack)
{  /*lint --e{715}*/
   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableKnapsack NULL
#endif




/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
static
DECL_LINCONSUPGD(linconsUpgdKnapsack)
{  /*lint --e{715}*/
   Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a knapsack constraint
    * - all variables must be binary
    * - all coefficients must be integral
    * - one of the sides must be infinite
    */
   upgrade = (nposbin + nnegbin == nvars)
      && (ncoeffspone + ncoeffsnone + ncoeffspint + ncoeffsnint == nvars)
      && (SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs));

   if( upgrade )
   {
      debugMessage("upgrading constraint <%s> to knapsack constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Knapsack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( SCIPcreateConsKnapsack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   }

   return SCIP_OKAY;
}
#endif




/*
 * constraint specific interface methods
 */

/** creates the handler for knapsack constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrKnapsack(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create knapsack constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeKnapsack, consInitKnapsack, consExitKnapsack, consSolstartKnapsack,
                  consDeleteKnapsack, consTransKnapsack, consInitlpKnapsack,
                  consSepaKnapsack, consEnfolpKnapsack, consEnfopsKnapsack, consCheckKnapsack, 
                  consPropKnapsack, consPresolKnapsack, consRescvarKnapsack,
                  consLockKnapsack, consUnlockKnapsack,
                  consActiveKnapsack, consDeactiveKnapsack, 
                  consEnableKnapsack, consDisableKnapsack,
                  conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdKnapsack, LINCONSUPGD_PRIORITY) );
#endif

   /* add knapsack constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a knapsack constraint */
RETCODE SCIPcreateConsKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             lhs,                /**< left hand side of constraint */
   Real             rhs,                /**< right hand side of constraint */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   errorMessage("method of knapsack constraint handler not implemented yet");
   abort(); /*lint --e{527} --e{715}*/

   /* find the knapsack constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("knapsack constraint handler not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                  local, modifiable, removeable) );

   return SCIP_OKAY;
}
