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

/**@file   cons_binpack.c
 * @brief  constraint handler for binpack constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_binpack.h"
#include "cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "binpack"
#define CONSHDLR_DESC          "bin packing constraints of the form  a^T x <= b*y, x, y binary"
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
DECL_CONSFREE(consFreeBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consFreeBinpack NULL
#endif


/** initialization method of constraint handler (called when problem solving starts) */
#if 0
static
DECL_CONSINIT(consInitBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consInitBinpack NULL
#endif


/** deinitialization method of constraint handler (called when problem solving exits) */
#if 0
static
DECL_CONSEXIT(consExitBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consExitBinpack NULL
#endif


/** frees specific constraint data */
#if 0
static
DECL_CONSDELETE(consDeleteBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consDeleteBinpack NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 0
static
DECL_CONSTRANS(consTransBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consTransBinpack NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
DECL_CONSINITLP(consInitlpBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consInitlpBinpack NULL
#endif


/** separation method of constraint handler */
#if 0
static
DECL_CONSSEPA(consSepaBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consSepaBinpack NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
DECL_CONSPROP(consPropBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consPropBinpack NULL
#endif


/** presolving method of constraint handler */
#if 0
static
DECL_CONSPRESOL(consPresolBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consPresolBinpack NULL
#endif


/** conflict variable resolving method of constraint handler */
#if 0
static
DECL_CONSRESCVAR(consRescvarBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consRescvarBinpack NULL
#endif


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
DECL_CONSACTIVE(consActiveBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consActiveBinpack NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
DECL_CONSDEACTIVE(consDeactiveBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consDeactiveBinpack NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
DECL_CONSENABLE(consEnableBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consEnableBinpack NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
DECL_CONSDISABLE(consDisableBinpack)
{
   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define consDisableBinpack NULL
#endif




/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
static
DECL_LINCONSUPGD(linconsUpgdBinpack)
{
   Bool upgrade;

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
      && (SCIPisInfinity(scip, -lhs) ^ SCIPisInfinity(scip, rhs)) )
   {
      Bool found;
      int i;

      found = FALSE;
      if( SCIPisInfinity(scip, -lhs) )
      {
         assert(!SCIPisInfinity(scip, rhs));

         for( i = 0; i < nvars && !found; ++i )
            found = SCIPisEQ(scip, rhs, ABS(vals[i]) + negcoeffsum);
      }
      else
      {
         assert(SCIPisInfinity(scip, rhs));

         for( i = 0; i < nvars && !found; ++i )
            found = SCIPisEQ(scip, lhs, -ABS(vals[i]) - negcoeffsum);
      }

      upgrade = found;
   }

   if( upgrade )
   {
      debugMessage("upgrading constraint <%s> to binpack constraint\n", SCIPconsGetName(cons));
      
      /* create the bin Binpack constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( SCIPcreateConsBinpack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, rhs,
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

/** creates the handler for binpack constraints and includes it in SCIP */
RETCODE SCIPincludeConsHdlrBinpack(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create binpack constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeBinpack, consInitBinpack, consExitBinpack,
                  consDeleteBinpack, consTransBinpack, consInitlpBinpack,
                  consSepaBinpack, consEnfolpBinpack, consEnfopsBinpack, consCheckBinpack, 
                  consPropBinpack, consPresolBinpack, consRescvarBinpack,
                  consLockBinpack, consUnlockBinpack,
                  consActiveBinpack, consDeactiveBinpack, 
                  consEnableBinpack, consDisableBinpack,
                  conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdBinpack, LINCONSUPGD_PRIORITY) );
#endif

   /* add binpack constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a binpack constraint */
RETCODE SCIPcreateConsBinpack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
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

   errorMessage("method of binpack constraint handler not implemented yet");
   abort();

   /* find the binpack constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("binpack constraint handler not found");
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
