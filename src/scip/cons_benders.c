/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_benders.c
 * @brief  constraint handler for benders decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/cons_benders.h"
#include "scip/cons_benderslp.h"
#include "scip/heur_trysol.h"
#include "scip/heuristics.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "benders"
#define CONSHDLR_DESC          "constraint handler to execute Benders' Decomposition"
#define CONSHDLR_ENFOPRIORITY        -1 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -5000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */


#define DEFAULT_CHECKEDSOLSSIZE      20 /**< initial size of the checked sols array */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   ncalls;             /**< the number of calls to the constraint handler. */
   int*                  checkedsols;        /**< an array of solutions that this constraint has already checked */
   int                   ncheckedsols;       /**< the number of checked solutions */
   int                   checkedsolssize;    /**< the size of the checked solutions array */
};

/*
 * Local methods
 */

/** constructs a new solution based upon the solutions to the Benders' decomposition subproblems */
static
SCIP_RETCODE constructValidSolution(
   SCIP*                 scip,               /**< the SCIP instance */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_SOL* newsol;
   SCIP_HEUR* heurtrysol;
   SCIP_BENDERS** benders;
   SCIP_VAR** auxiliaryvars;
   int nactivebenders;
   int nsubproblems;
   int i;
   int j;
   SCIP_Bool success;


   /* don't propose new solutions if not in presolve or solving */
   if( SCIPgetStage(scip) < SCIP_STAGE_INITPRESOLVE || SCIPgetStage(scip) >= SCIP_STAGE_SOLVED )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   benders = SCIPgetBenders(scip);
   nactivebenders = SCIPgetNActiveBenders(scip);

   /* if the solution is NULL, then we create the solution from the LP sol */
   if( sol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );

   /* checking the size of the checkedsols array and extending it is there is not enough memory */
   assert(conshdlrdata->ncheckedsols <= conshdlrdata->checkedsolssize);
   if( conshdlrdata->ncheckedsols + 1 > conshdlrdata->checkedsolssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->ncheckedsols + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->checkedsols, conshdlrdata->checkedsolssize, newsize) );
      conshdlrdata->checkedsolssize = newsize;
   }
   assert(conshdlrdata->ncheckedsols + 1 <= conshdlrdata->checkedsolssize);

   /* recording the solution number to avoid checking the solution again */
   conshdlrdata->checkedsols[conshdlrdata->ncheckedsols] = SCIPsolGetIndex(newsol);
   conshdlrdata->ncheckedsols++;

   /* looping through all Benders' decompositions to construct the new solution */
   for( i = 0; i < nactivebenders; i++ )
   {
      /* getting the auxiliary variables and the number of subproblems from the Benders' decomposition structure */
      auxiliaryvars = SCIPbendersGetAuxiliaryVars(benders[i]);
      nsubproblems = SCIPbendersGetNSubproblems(benders[i]);

      /* setting the auxiliary variable in the new solution */
      for( j = 0; j < nsubproblems; j++ )
         SCIP_CALL( SCIPsetSolVal(scip, newsol, auxiliaryvars[j], SCIPbendersGetSubprobObjval(benders[i], j)) );
   }

   /* getting the try solution heuristic */
   heurtrysol = SCIPfindHeur(scip, "trysol");

   /* passing the new solution to the trysol heuristic  */
   SCIP_CALL( SCIPcheckSol(scip, newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
   if ( success )
   {
      SCIP_CALL( SCIPheurPassSolAddSol(scip, heurtrysol, newsol) );
      SCIPdebugMsg(scip, "Creating solution was successful.\n");
   }
#ifdef SCIP_DEBUG
   else
   {
      /* the solution might not be feasible, because of additional constraints */
      SCIPdebugMsg(scip, "Creating solution was not successful.\n");
   }
#endif

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** the methods for the enforcement of solutions */
SCIP_RETCODE SCIPconsBendersEnforceSolutions(
   SCIP*                 scip,               /**< the SCIP instance */
   SCIP_SOL*             sol,                /**< the primal solution to enforce, or NULL for the current LP/pseudo sol */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_RESULT*          result,             /**< the result of the enforcement */
   SCIP_BENDERSENFOTYPE  type,               /**< the type of solution being enforced */
   SCIP_Bool             checkint            /**< should the integer solution be checked by the subproblems */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_BENDERS** benders;
   SCIP_Bool infeasible;
   SCIP_Bool auxviol;
   int nactivebenders;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;
   infeasible = FALSE;
   auxviol = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   benders = SCIPgetBenders(scip);
   nactivebenders = SCIPgetNActiveBenders(scip);

   for( i = 0; i < nactivebenders; i++ )
   {
      switch( type )
      {
         case SCIP_BENDERSENFOTYPE_LP:
            if( SCIPbendersCutLP(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], NULL, result, &infeasible, &auxviol, type, checkint) );
            }
            break;
         case SCIP_BENDERSENFOTYPE_RELAX:
            if( SCIPbendersCutRelaxation(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], sol, result, &infeasible, &auxviol, type, checkint) );
            }
            break;
         case SCIP_BENDERSENFOTYPE_PSEUDO:
            if( SCIPbendersCutPseudo(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], NULL, result, &infeasible, &auxviol, type, checkint) );
            }
            break;
         case SCIP_BENDERSENFOTYPE_CHECK:
            SCIPwarningMessage(scip, "The conscheck callback is not supported\n");
            break;
         default:
            break;
      }

      /* The decompositions are checked until one is found that is not feasible. Not being feasible could mean that
       * infeasibility of the original problem has been proven or a constraint has been added. If the result
       * SCIP_DIDNOTRUN is returned, then the next decomposition is checked */
      if( (*result) != SCIP_FEASIBLE && (*result) != SCIP_DIDNOTRUN )
         break;
   }

   /* if the constraint handler was called with an integer feasible solution, then a feasible solution can be proposed */
   if( checkint )
   {
      /* in the case that the problem is feasible, this means that all subproblems are feasible. The auxiliary variables
       * still need to be updated. This is done by constructing a valid solution. */
      if( (*result) == SCIP_FEASIBLE && auxviol )
      {
         if( type == SCIP_BENDERSENFOTYPE_PSEUDO )
         {
            if( !SCIPsolIsOriginal(sol) )
               SCIP_CALL( constructValidSolution(scip, conshdlr, sol) );
         }

         (*result) = SCIP_INFEASIBLE;
      }
   }

   /* if no Benders' decomposition were run, then the result is returned as SCIP_FEASIBLE. The SCIP_DIDNOTRUN result
    * indicates that no subproblems were checked */
   if( (*result) == SCIP_DIDNOTRUN )
      (*result) = SCIP_FEASIBLE;

   conshdlrdata->ncalls++;

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBenders)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeConshdlrBenders(scip, FALSE) );

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* freeing the constraint handler data */
   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   conshdlrdata->checkedsolssize = DEFAULT_CHECKEDSOLSSIZE;
   conshdlrdata->ncheckedsols = 0;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->checkedsols, conshdlrdata->checkedsolssize) );

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* freeing the checked sols array */
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->checkedsols, conshdlrdata->checkedsolssize);

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBenders)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, NULL, conshdlr, result, SCIP_BENDERSENFOTYPE_LP, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBenders)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, sol, conshdlr, result, SCIP_BENDERSENFOTYPE_RELAX, TRUE) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBenders)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, NULL, conshdlr, result, SCIP_BENDERSENFOTYPE_PSEUDO, TRUE) );

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
/*  This function checks the feasibility of the Benders' decomposition master problem. In the case that the problem is
 *  feasible, then the auxiliary variables must be updated with the subproblem objective function values. It is not
 *  possible to simply update the auxiliary variable values, so a new solution is created. */
static
SCIP_DECL_CONSCHECK(consCheckBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_BENDERS** benders;
   int nactivebenders;
   int solindex;
   int i;
   SCIP_Bool performcheck;
   SCIP_Bool infeasible;
   SCIP_Bool auxviol;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;
   performcheck = TRUE;
   infeasible = FALSE;
   auxviol = FALSE;


   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   benders = SCIPgetBenders(scip);
   nactivebenders = SCIPgetNActiveBenders(scip);

   /* checking if the solution was constructed by this constraint handler */
   solindex = SCIPsolGetIndex(sol);
   for( i = 0; i < conshdlrdata->ncheckedsols; i++ )
   {
      if( conshdlrdata->checkedsols[i] == solindex )
      {
         conshdlrdata->checkedsols[0] = conshdlrdata->checkedsols[conshdlrdata->ncheckedsols - 1];
         conshdlrdata->ncheckedsols--;

         performcheck = FALSE;
         break;
      }
   }

   /* if the solution has not been checked before, then we must perform the check */
   if( performcheck && nactivebenders > 0 )
   {
      for( i = 0; i < nactivebenders; i++ )
      {
         SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], sol, result, &infeasible, &auxviol,
               SCIP_BENDERSENFOTYPE_CHECK, TRUE) );

         /* in the case of multiple Benders' decompositions, the subproblems are solved until a constriant is added or
          * infeasibility is proven. So if the result is not SCIP_FEASIBLE, then the loop is exited */
         if( (*result) != SCIP_FEASIBLE )
            break;
      }

      /* in the case that the problem is feasible, this means that all subproblems are feasible. The auxiliary variables
       * still need to be updated. This is done by constructing a valid solution. */
      if( (*result) == SCIP_FEASIBLE )
      {
         if( auxviol )
         {
            if( !SCIPsolIsOriginal(sol) )
               SCIP_CALL( constructValidSolution(scip, conshdlr, sol) );

            if( printreason )
               SCIPmessagePrintInfo(SCIPgetMessagehdlr(scip), "all subproblems are feasible but there is a violation in the auxiliary variables\n");

            (*result) = SCIP_INFEASIBLE;
         }
      }

      /* if no Benders' decomposition were run, then the result is returned as SCIP_FEASIBLE. The SCIP_DIDNOTRUN result
       * indicates that no subproblems were checked */
      if( (*result) == SCIP_DIDNOTRUN )
         (*result) = SCIP_FEASIBLE;
   }

   conshdlrdata->ncalls++;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBenders)
{  /*lint --e{715}*/
   //SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   //SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for benders constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             twophase            /**< should the two phase method be used? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create benders constraint handler data */
   conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->ncalls = 0;

   conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBenders, consEnfopsBenders, consCheckBenders, consLockBenders,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitBenders) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitBenders) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBenders, NULL) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBenders) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBenders) );

   if( twophase )
      SCIP_CALL( SCIPincludeConshdlrBenderslp(scip) );

   return SCIP_OKAY;
}

/** creates and captures a benders constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBenders(
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsBenders() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the benders constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("benders constraint handler not found\n");
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

/** creates and captures a benders constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicBenders(
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
   SCIP_CALL( SCIPcreateConsBenders(scip, cons, name, nvars, vars, coefs, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
