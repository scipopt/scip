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

/**@file   cons_eqknapsack.c
 * @brief  constraint handler for equality knapsack constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_eqknapsack.h"
#include "cons_linear.h"


#define CONSHDLR_NAME          "eqknapsack"
#define CONSHDLR_DESC          "equality knapsack constraints"
#define CONSHDLR_SEPAPRIORITY   +400100
#define CONSHDLR_ENFOPRIORITY   +400100
#define CONSHDLR_CHECKPRIORITY  -400100
#define CONSHDLR_SEPAFREQ            -1
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS         TRUE

#define LINCONSUPGD_PRIORITY    +400100




/*
 * Local methods
 */




/*
 * Callback methods of constraint handler
 */




/*
 * Presolving
 */




/*
 * Linear constraint upgrading
 */

static
DECL_LINCONSUPGD(linconsUpgdEqknapsack)
{
   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a equality knapsack constraint
    * - all variables must be binary
    * - all coefficients must be integral
    * - both sides must be equal
    */
   if( nposbin + nnegbin == nvars
      && ncoeffspone + ncoeffsnone + ncoeffspint + ncoeffsnint == nvars
      && SCIPisEQ(scip, lhs, rhs) )
   {
      debugMessage("upgrading constraint <%s> to equality knapsack constraint\n", SCIPconsGetName(cons));
      
      /* create the equality knapsack constraint (an automatically upgraded constraint is always unmodifiable) */
      CHECK_OKAY( SCIPcreateConsEqknapsack(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, rhs,
                     SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
                     SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
                     local, FALSE, removeable) );
   }

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for equality knapsack constraints and includes it in SCIP */
RETCODE SCIPincludeConsHdlrEqknapsack(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConsHdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ,
                  CONSHDLR_NEEDSCONS,
                  NULL, NULL, NULL,
                  NULL, NULL, 
                  NULL, NULL, NULL, NULL, NULL, NULL,
                  NULL, NULL,
                  NULL) );

   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdEqknapsack, LINCONSUPGD_PRIORITY) );

   return SCIP_OKAY;
}

/** creates and captures a equality knapsack constraint */
RETCODE SCIPcreateConsEqknapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR**            vars,               /**< array with variables of constraint entries */
   Real*            vals,               /**< array with coefficients of constraint entries */
   Real             rhs,                /**< right hand side of row */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   )
{
   CONSHDLR* conshdlr;
   CONSDATA* consdata;

   errorMessage("equality knapsack constraint handler not implemented yet");
   abort();

   /* find the linear constraint handler */
   conshdlr = SCIPfindConsHdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("equality knapsack constraint handler not found");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;

   /* create constraint */
   CHECK_OKAY( SCIPcreateCons(scip, cons, name, conshdlr, consdata, separate, enforce, check, propagate) );

   return SCIP_OKAY;
}
