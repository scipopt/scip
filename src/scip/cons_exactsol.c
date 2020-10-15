/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_exactsol.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for ensuring that primal solution is exact
 * @author Antonia Chmiela
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_exactsol.h"
#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/pub_var.h"
#include "scip/rational.h"
//#include "scip/struct_rational.h"
#include "scip/cons.h"
#include "scip/scip_exact.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_sol.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/struct_scip.h"
#include "scip/tree.h"



/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "exactsol"
#define CONSHDLR_DESC          "constraint to ensure that primal solutions report back exact solutions"
#define CONSHDLR_ENFOPRIORITY  -9999999 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
//#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
//#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
//#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

//#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
//#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
//#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

//#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
//#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */




/* TODO: (optional) enable linear or nonlinear constraint upgrading */
#if 0
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#define LINCONSUPGD_PRIORITY          0 /**< priority of the constraint handler for upgrading of linear constraints */
#define NONLINCONSUPGD_PRIORITY       0 /**< priority of the constraint handler for upgrading of nonlinear constraints */
#endif


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for ExactSol constraints */
//struct SCIP_ConsData
//{
//};

/** constraint handler data */
//struct SCIP_ConshdlrData
//{
//};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of constraint handler
 */


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExactSol)
{  /*lint --e{715}*/

   /* returning feasible since we can't enforce anything */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxExactSol)
{  /*lint --e{715}*/

   /* returning feasible since we can't enforce anything */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExactSol)
{  /*lint --e{715}*/

   /* returning feasible since we can't enforce anything */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** TODO */
/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExactSol)
{  /*lint --e{715}*/

   SCIP_VAR** vars;
   SCIP_Bool lperror;
   int nintvars;
   int nvars;
   int i;
   int c;
#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* if the solution doesn't come from a heuristic, ignore it */
   if( SCIPsolGetType(sol) != SCIP_SOLTYPE_HEUR )
      return SCIP_OKAY;

   /* if we are not solving exactly, we have nothing to check */
   if( !SCIPisExactSolve(scip) )
      return SCIP_OKAY;

   /** if we're already in exact diving mode, we already computed an exact solution with this constraint handler and
     * are checking if it's actually feasible
     */
   if( SCIPlpExactDiving(scip->lpexact) )
      return SCIP_OKAY;

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   if( !SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIP_OKAY;

   if( !SCIPtreeHasCurrentNodeLP(scip->tree))
      return SCIP_OKAY;

   if( !SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_FOCUSNODE )
      return SCIP_OKAY;

   if( SCIPisExactSol(scip, sol) )
      return SCIP_OKAY;

   /* check if solution is floating point feasible */
   for( c = 0; c < nconss && *result == SCIP_FEASIBLE; ++c )
   {
      SCIP_Real activity;
      SCIP_ROW* row;

      /* get row corresponding to constraint */
      row = SCIPconsGetRow(scip, conss[c]);
      assert(row != NULL);

      /* get row activity */
      activity = SCIPgetRowSolActivity(scip, row, sol);

      /* check if the constraint is violated */
      if( SCIPisFeasLT(scip, activity, SCIProwGetLhs(row)) || SCIPisFeasGT(scip, activity, SCIProwGetRhs(row)) )
         *result = SCIP_INFEASIBLE;
   }

   if( *result == SCIP_INFEASIBLE )
      return SCIP_OKAY;

   *result = SCIP_INFEASIBLE;

   /* start exact diving */
   SCIP_CALL( SCIPstartExactDive(scip) );

   /* set the bounds of the variables: fixed for integers, global bounds for continuous */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

   for( i = 0; i < nintvars; ++i )
   {
      if( SCIPvarGetStatusExact(vars[i]) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_Real solval;
         solval = SCIPgetSolVal(scip, sol, vars[i]);
         //printf(" (exactsol) lb = %f, ub = %f of var %s \n", RatApproxReal(SCIPvarGetLbLocalExact(vars[i])), RatApproxReal(SCIPvarGetUbLocalExact(vars[i])), SCIPvarGetName(vars[i]));
         assert(RatIsLE(SCIPvarGetLbLocalExact(vars[i]), SCIPvarGetUbLocalExact(vars[i])));

         /* check if solution value is supposed to be integral */
         if( SCIPisIntegral(scip, solval) )
         {
            SCIP_Rational* newbound;

            RatCreateBuffer(SCIPbuffer(scip), &newbound);

            /* create rational solval and round it to the nearest integer */
            RatSetReal(newbound, solval);
            RatRound(newbound, newbound, SCIP_ROUND_NEAREST);

            //printf("solval = %f \n",solval );
            //printf("fp = <%f>, exact = <%f> ---> %d \n", SCIPsetRound(scip->set, solval), RatApproxReal(newbound), RatIsApproxEqualReal(scip->set, newbound, SCIPsetRound(scip->set, solval)));

            SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], SCIPsetRound(scip->set, solval)) );
            SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], SCIPsetRound(scip->set, solval)) );

            SCIP_CALL( SCIPchgVarLbExactDive(scip, vars[i], newbound) );
            SCIP_CALL( SCIPchgVarUbExactDive(scip, vars[i], newbound) );

            RatFreeBuffer(SCIPbuffer(scip), &newbound);
         }
      }
   }

   /* solve LP */

   /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a constraint
    * handler. Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP
    * will stop.
    */
#ifdef NDEBUG
   retstat = SCIPsolveExactDiveLP(scip, -1, &lperror, NULL);
   if( retstat != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while solving LP in Exactsol Constraint Handler; Exact LP solve terminated with code <%d>\n",retstat);
   }
#else
   SCIP_CALL( SCIPsolveExactDiveLP(scip, -1, &lperror, NULL) );
#endif

   /* check if this is a feasible solution */
   if( !lperror && SCIPgetLPExactSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_SOL* exactsol;
      SCIP_Bool foundsol;

      SCIP_CALL( SCIPsolCreateLPSolExact(&exactsol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
          scip->tree, scip->lpexact, NULL) );
      SCIPsolSetHeur(exactsol, SCIPsolGetHeur(sol));
      SCIP_CALL( SCIPtrySolFreeExact(scip, &exactsol, FALSE, FALSE, FALSE, FALSE, TRUE, &foundsol) );

      /* check solution for feasibility, and add it to solution store if possible
       * neither integrality nor feasibility of LP rows has to be checked, because this is already
       * done in the intshifting heuristic itself and due to the LP resolve
       */
      //SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, FALSE, &stored) );

      if( foundsol )
      {
         *result = SCIP_FEASIBLE;
      }
   }

   /* terminate exact diving */
   SCIP_CALL( SCIPendExactDive(scip) );

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExactSol)
{  /*lint --e{715}*/

   /* do nothing since we are not handling constraints */
   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for ExactSol constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExactSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create ExactSol constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

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
         conshdlrCopyExactSol,
         consFreeExactSol, consInitExactSol, consExitExactSol,
         consInitpreExactSol, consExitpreExactSol, consInitsolExactSol, consExitsolExactSol,
         consDeleteExactSol, consTransExactSol, consInitlpExactSol,
         consSepalpExactSol, consSepasolExactSol, consEnfolpExactSol, consEnforelaxExactSol, consEnfopsExactSol, consCheckExactSol,
         consPropExactSol, consPresolExactSol, consRespropExactSol, consLockExactSol,
         consActiveExactSol, consDeactiveExactSol,
         consEnableExactSol, consDisableExactSol, consDelvarsExactSol,
         consPrintExactSol, consCopyExactSol, consParseExactSol,
         consGetVarsExactSol, consGetNVarsExactSol, consGetDiveBdChgsExactSol, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpExactSol, consEnfopsExactSol, consCheckExactSol, consLockExactSol,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* add ExactSol constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */
#endif
   return SCIP_OKAY;
}

/** creates and captures a ExactSol constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExactSol(
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsExactSol() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of ExactSol constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the ExactSol constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("ExactSol constraint handler not found\n");
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

/** creates and captures a ExactSol constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicExactSol(
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
   SCIP_CALL( SCIPcreateConsExactSol(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
