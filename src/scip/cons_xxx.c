/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_xxx.c,v 1.16 2004/05/03 08:13:09 bzfpfend Exp $"

/**@file   cons_xxx.c
 * @brief  constraint handler for xxx constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cons_xxx.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "xxx"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_SEPAPRIORITY         0
#define CONSHDLR_ENFOPRIORITY         0
#define CONSHDLR_CHECKPRIORITY        0
#define CONSHDLR_SEPAFREQ            -1
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE

/* TODO: (optional) enable linear constraint upgrading */
#if 0
#include "cons_linear.h"

#define LINCONSUPGD_PRIORITY          0
#endif




/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for xxx constraints */
struct ConsData
{
};

/** constraint handler data */
struct ConshdlrData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
DECL_CONSFREE(consFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeXxx NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
DECL_CONSINIT(consInitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitXxx NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
DECL_CONSEXIT(consExitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitXxx NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
DECL_CONSINITPRE(consInitpreXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreXxx NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
DECL_CONSEXITPRE(consExitpreXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreXxx NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
DECL_CONSINITSOL(consInitsolXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolXxx NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
DECL_CONSEXITSOL(consExitsolXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolXxx NULL
#endif


/** frees specific constraint data */
#if 0
static
DECL_CONSDELETE(consDeleteXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteXxx NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 0
static
DECL_CONSTRANS(consTransXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransXxx NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
DECL_CONSINITLP(consInitlpXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpXxx NULL
#endif


/** separation method of constraint handler */
#if 0
static
DECL_CONSSEPA(consSepaXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepaXxx NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
DECL_CONSPROP(consPropXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropXxx NULL
#endif


/** presolving method of constraint handler */
#if 0
static
DECL_CONSPRESOL(consPresolXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolXxx NULL
#endif


/** conflict variable resolving method of constraint handler */
#if 0
static
DECL_CONSRESCVAR(consRescvarXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRescvarXxx NULL
#endif


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
DECL_CONSACTIVE(consActiveXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveXxx NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
DECL_CONSDEACTIVE(consDeactiveXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveXxx NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
DECL_CONSENABLE(consEnableXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableXxx NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
DECL_CONSDISABLE(consDisableXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableXxx NULL
#endif




/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a xxx constraint */
static
DECL_LINCONSUPGD(linconsUpgdXxx)
{  /*lint --e{715}*/
   Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to xxx constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      debugMessage("upgrading constraint <%s> to xxx constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Xxx constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( SCIPcreateConsXxx(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
                     SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemoveable(cons)) );
   }

   return SCIP_OKAY;
}
#endif




/*
 * constraint specific interface methods
 */

/** creates the handler for xxx constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create xxx constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeXxx, consInitXxx, consExitXxx, 
                  consInitpreXxx, consExitpreXxx, consInitsolXxx, consExitsolXxx,
                  consDeleteXxx, consTransXxx, consInitlpXxx,
                  consSepaXxx, consEnfolpXxx, consEnfopsXxx, consCheckXxx, 
                  consPropXxx, consPresolXxx, consRescvarXxx,
                  consLockXxx, consUnlockXxx,
                  consActiveXxx, consDeactiveXxx, 
                  consEnableXxx, consDisableXxx,
                  conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdXxx, LINCONSUPGD_PRIORITY) );
#endif

   /* add xxx constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a xxx constraint */
RETCODE SCIPcreateConsXxx(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of variables in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            coefs,              /**< array with coefficients of constraint entries */
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsXxx() call, if you don't need all the information */

   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   errorMessage("method of xxx constraint handler not implemented yet\n");
   abort(); /*lint --e{527} --e{715}*/

   /* find the xxx constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("xxx constraint handler not found\n");
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
