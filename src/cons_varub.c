/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_varub.c,v 1.10 2003/11/21 10:35:34 bzfpfend Exp $"

/**@file   cons_varub.c
 * @brief  constraint handler for varub constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_varub.h"
#include "cons_linear.h"


/* constraint handler properties */
#define CONSHDLR_NAME          "varub"
#define CONSHDLR_DESC          "variable upper bounds of the form  y <= a*x, x binary"
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
DECL_CONSFREE(consFreeVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeVarub NULL
#endif


/** initialization method of constraint handler (called when problem solving starts) */
#if 0
static
DECL_CONSINIT(consInitVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitVarub NULL
#endif


/** deinitialization method of constraint handler (called when problem solving exits) */
#if 0
static
DECL_CONSEXIT(consExitVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitVarub NULL
#endif


/** solving start notification method of constraint handler (called when presolving was finished) */
#if 0
static
DECL_CONSSOLSTART(consSolstartVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSolstartVarub NULL
#endif


/** frees specific constraint data */
#if 0
static
DECL_CONSDELETE(consDeleteVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteVarub NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */ 
#if 0
static
DECL_CONSTRANS(consTransVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransVarub NULL
#endif


/** LP initialization method of constraint handler */
#if 0
static
DECL_CONSINITLP(consInitlpVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpVarub NULL
#endif


/** separation method of constraint handler */
#if 0
static
DECL_CONSSEPA(consSepaVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepaVarub NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
DECL_CONSENFOPS(consEnfopsVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
DECL_CONSPROP(consPropVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropVarub NULL
#endif


/** presolving method of constraint handler */
#if 0
static
DECL_CONSPRESOL(consPresolVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolVarub NULL
#endif


/** conflict variable resolving method of constraint handler */
#if 0
static
DECL_CONSRESCVAR(consRescvarVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRescvarVarub NULL
#endif


/** variable rounding lock method of constraint handler */
static
DECL_CONSLOCK(consLockVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** variable rounding unlock method of constraint handler */
static
DECL_CONSUNLOCK(consUnlockVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
DECL_CONSACTIVE(consActiveVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveVarub NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
DECL_CONSDEACTIVE(consDeactiveVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveVarub NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
DECL_CONSENABLE(consEnableVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableVarub NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
DECL_CONSDISABLE(consDisableVarub)
{  /*lint --e{715}*/
   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableVarub NULL
#endif




/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
static
DECL_LINCONSUPGD(linconsUpgdVarub)
{  /*lint --e{715}*/
   Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a variable upper bound constraint
    * - there is one binary variable and one non-binary variable
    * - if the non-binary variable coefficient is positive, the right hand side is zero and the left hand side is infinite,
    *   if the non-binary variable coefficient is negative, the left hand side is zero and the right hand side is infinite
    */
   upgrade = (nvars == 2)
      && (nposbin + nnegbin == 1 && nposint + nnegint + nposimpl + nnegimpl + nposcont + nnegcont == 1)
      && ((nposint + nposimpl + nposcont == 1 && SCIPisZero(scip, rhs) && SCIPisInfinity(scip, -lhs))
         || (nnegint + nnegimpl + nnegcont == 1 && SCIPisZero(scip, lhs) && SCIPisInfinity(scip, rhs)));

   if( upgrade )
   {
      VAR* var;
      VAR* switchvar;
      Real val;

      debugMessage("upgrading constraint <%s> to varub constraint\n", SCIPconsGetName(cons));

      /* find out the variable with variable bound, and the bound switching variable */
      if( SCIPvarGetType(vars[0]) == SCIP_VARTYPE_BINARY )
      {
         switchvar = vars[0];
         var = vars[1];
         val = -vals[0]/vals[1];
      }
      else
      {
         switchvar = vars[1];
         var = vars[0];
         val = -vals[1]/vals[0];
      }

      /* create the bin Varub constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      CHECK_OKAY( SCIPcreateConsVarub(scip, upgdcons, SCIPconsGetName(cons), nvars, var, switchvar, val,
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

/** creates the handler for varub constraints and includes it in SCIP */
RETCODE SCIPincludeConshdlrVarub(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create varub constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeVarub, consInitVarub, consExitVarub, consSolstartVarub,
                  consDeleteVarub, consTransVarub, consInitlpVarub,
                  consSepaVarub, consEnfolpVarub, consEnfopsVarub, consCheckVarub, 
                  consPropVarub, consPresolVarub, consRescvarVarub,
                  consLockVarub, consUnlockVarub,
                  consActiveVarub, consDeactiveVarub, 
                  consEnableVarub, consDisableVarub,
                  conshdlrdata) );

#ifdef LINCONSUPGD_PRIORITY
   /* include the linear constraint upgrade in the linear constraint handler */
   CHECK_OKAY( SCIPincludeLinconsUpgrade(scip, linconsUpgdVarub, LINCONSUPGD_PRIORITY) );
#endif

   /* add varub constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a varub constraint */
RETCODE SCIPcreateConsVarub(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              len,                /**< number of nonzeros in the constraint */
   VAR*             var,                /**< variable that has variable bound */
   VAR*             switchvar,          /**< binary variable to activate bound */
   Real             val,                /**< bound value */
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

   errorMessage("method of varub constraint handler not implemented yet\n");
   abort(); /*lint --e{527} --e{715}*/

   /* find the varub constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      errorMessage("varub constraint handler not found\n");
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
