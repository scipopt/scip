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
#pragma ident "@(#) $Id: cons_integral.c,v 1.24 2004/04/27 15:49:58 bzfpfend Exp $"

/**@file   cons_integral.c
 * @brief  constraint handler for the integrality constraint
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "cons_integral.h"


#define CONSHDLR_NAME          "integral"
#define CONSHDLR_DESC          "integrality constraint"
#define CONSHDLR_SEPAPRIORITY         0
#define CONSHDLR_ENFOPRIORITY         0
#define CONSHDLR_CHECKPRIORITY        0
#define CONSHDLR_SEPAFREQ            -1
#define CONSHDLR_PROPFREQ            -1
#define CONSHDLR_NEEDSCONS        FALSE /**< the constraint handler is called without constraints */



/*
 * Callback methods
 */

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#define consFreeIntegral NULL


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitIntegral NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitIntegral NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolIntegral NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolIntegral NULL


/** frees specific constraint data */
#define consDeleteIntegral NULL


/** transforms constraint data into data belonging to the transformed problem */ 
#define consTransIntegral NULL


/** LP initialization method of constraint handler */
#define consInitlpIntegral NULL


/** separation method of constraint handler */
#define consSepaIntegral NULL


/** constraint enforcing method of constraint handler for LP solutions */
static
DECL_CONSENFOLP(consEnfolpIntegral)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);

   debugMessage("Enfolp method of integrality constraint\n");

   /* call branching methods */
   CHECK_OKAY( SCIPbranchLP(scip, result) );

   /* if no branching was done, the LP solution was not fractional */
   if( *result == SCIP_DIDNOTRUN )
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
#define consEnfopsIntegral NULL


/** feasibility check method of constraint handler for integral solutions */
static
DECL_CONSCHECK(consCheckIntegral)
{  /*lint --e{715}*/
   VAR** vars;
   Real solval;
   int nbin;
   int nint;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   debugMessage("Check method of integrality constraint\n");

   CHECK_OKAY( SCIPgetVarsData(scip, &vars, NULL, &nbin, &nint, NULL, NULL) );

   *result = SCIP_FEASIBLE;

   if( checkintegrality )
   {
      for( v = 0; v < nbin + nint && *result == SCIP_FEASIBLE; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, vars[v]);
         if( !SCIPisIntegral(scip, solval) )
            *result = SCIP_INFEASIBLE;
      }
   }
#ifndef NDEBUG
   else
   {
      for( v = 0; v < nbin + nint; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, vars[v]);
         assert(SCIPisIntegral(scip, solval));
      }
   }
#endif

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#define consPropIntegral NULL


/** presolving method of constraint handler */
#define consPresolIntegral NULL


/** conflict variable resolving method of constraint handler */
#define consRescvarIntegral NULL


/** variable rounding lock method of constraint handler */
#define consLockIntegral NULL


/** variable rounding unlock method of constraint handler */
#define consUnlockIntegral NULL


/** constraint activation notification method of constraint handler */
#define consActiveIntegral NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveIntegral NULL


/** constraint enabling notification method of constraint handler */
#define consEnableIntegral NULL


/** constraint disabling notification method of constraint handler */
#define consDisableIntegral NULL




/*
 * constraint specific interface methods
 */

/** creates the handler for integrality constraint and includes it in SCIP */
RETCODE SCIPincludeConshdlrIntegral(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   CONSHDLRDATA* conshdlrdata;

   /* create integral constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   CHECK_OKAY( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
                  CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
                  CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_NEEDSCONS,
                  consFreeIntegral, consInitIntegral, consExitIntegral, consInitsolIntegral, consExitsolIntegral,
                  consDeleteIntegral, consTransIntegral, consInitlpIntegral,
                  consSepaIntegral, consEnfolpIntegral, consEnfopsIntegral, consCheckIntegral, 
                  consPropIntegral, consPresolIntegral, consRescvarIntegral,
                  consLockIntegral, consUnlockIntegral,
                  consActiveIntegral, consDeactiveIntegral, 
                  consEnableIntegral, consDisableIntegral,
                  conshdlrdata) );

   return SCIP_OKAY;
}
