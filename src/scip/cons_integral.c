/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_integral.c
 * @brief  constraint handler for the integrality constraint
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>

#include "scip/cons_integral.h"


#define CONSHDLR_NAME          "integral"
#define CONSHDLR_DESC          "integrality constraint"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/*
 * Callback methods
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyIntegral)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrIntegral(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

#define consCopyIntegral NULL

#define consEnfopsIntegral NULL

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpIntegral)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);

   SCIPdebugMessage("Enfolp method of integrality constraint: %d fractional variables\n", SCIPgetNLPBranchCands(scip));

   /* if the root LP is unbounded, we want to terminate with UNBOUNDED or INFORUNBOUNDED,
    * depending on whether we are able to construct an integral solution; in any case we do not want to branch
    */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      if( SCIPgetNLPBranchCands(scip) == 0 )
         *result = SCIP_FEASIBLE;
      else
         *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   /* call branching methods */
   SCIP_CALL( SCIPbranchLP(scip, result) );

   /* if no branching was done, the LP solution was not fractional */
   if( *result == SCIP_DIDNOTRUN )
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckIntegral)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real solval;
   int nallinteger;
   int ninteger;
   int nbin;
   int nint;
   int nimpl;
   int v;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);

   SCIPdebugMessage("Check method of integrality constraint (checkintegrality=%u)\n", checkintegrality);

   SCIP_CALL( SCIPgetSolVarsData(scip, sol, &vars, NULL, &nbin, &nint, &nimpl, NULL) );

   *result = SCIP_FEASIBLE;

   ninteger = nbin + nint;

   if( checkintegrality )
   {
      for( v = 0; v < ninteger; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, vars[v]);
         if( !SCIPisFeasIntegral(scip, solval) )
         {
            *result = SCIP_INFEASIBLE;

            if( printreason )
            {
               SCIPinfoMessage(scip, NULL, "violation: integrality condition of variable <%s> = %.15g\n", 
                  SCIPvarGetName(vars[v]), solval);
            }
            break;
         }
      }
   }
#ifndef NDEBUG
   else
   {
      for( v = 0; v < ninteger; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, vars[v]);
         assert(SCIPisFeasIntegral(scip, solval));
      }
   }
#endif

   nallinteger = ninteger + nimpl;
   for( v = ninteger; v < nallinteger; ++v )
   {
      solval = SCIPgetSolVal(scip, sol, vars[v]);
      if( !SCIPisFeasIntegral(scip, solval) )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIPinfoMessage(scip, NULL, "violation: integrality condition of implicit integral variable <%s> = %.15g\n",
               SCIPvarGetName(vars[v]), solval);
         }
         break;
      }
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockIntegral)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for integrality constraint and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create integral constraint handler data */
   conshdlrdata = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpIntegral, consEnfopsIntegral, consCheckIntegral, consLockIntegral,
         conshdlrdata) );

   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyIntegral, consCopyIntegral) );

   return SCIP_OKAY;
}
