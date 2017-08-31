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


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "benders"
#define CONSHDLR_DESC          "constraint handler to execute Benders' Decomposition"
#define CONSHDLR_ENFOPRIORITY        -1 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -5000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */


#define DEFAULT_CHECKEDSOLSSIZE     20    /**< the initial size of the checked sols array */

/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for benders constraints */
//struct SCIP_ConsData
//{
//};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   ncalls;             /**< the number of calls to the constraint handler. */
   int*                  checkedsols;        /**< an array of solutions that this constraint has already checked */
   int                   ncheckedsols;       /**< the number of checked solutions */
   int                   checkedsolssize;    /**< the size of the checked solutions array */
};

/** the methods for the enforcement of solutions */
SCIP_RETCODE SCIPconsBendersEnforceSolutions(
   SCIP*                 scip,               /**< the SCIP instance */
   SCIP_SOL*             sol,                /**< the primal solution to enforce, or NULL for the current LP/pseudo sol */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_RESULT*          result,             /**< the result of the enforcement */
   SCIP_BENDERSENFOTYPE  type                /**< the type of solution being enforced */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_BENDERS** benders;
   SCIP_Bool infeasible;
   int nactivebenders;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;
   infeasible = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   benders = SCIPgetBenders(scip);
   nactivebenders = SCIPgetNActiveBenders(scip);

   for( i = 0; i < nactivebenders; i++ )
   {
      switch( type )
      {
         case LP:
            if( SCIPbendersCutLP(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], NULL, result, &infeasible, type) );
            }
            break;
         case RELAX:
            if( SCIPbendersCutRelaxation(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], sol, result, &infeasible, type) );
            }
            break;
         case PSEUDO:
            if( SCIPbendersCutPseudo(benders[i]) )
            {
               SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], NULL, result, &infeasible, type) );
            }
            break;
         case CHECK:
            SCIPwarningMessage(scip, "The conscheck callback is not supported\n");
         default:
            break;
      }
   }

   conshdlrdata->ncalls++;

   return SCIP_OKAY;
}

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

   /* if the solution is NULL, then we create the solution from the LP sol. */
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

/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 1
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBenders)
{  /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeConshdlrBenders(scip, FALSE) );

   return SCIP_OKAY;
}
#else
#define conshdlrCopyBenders NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeBenders)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* freeing the constrains handler data */
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


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreBenders NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreBenders NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolBenders NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolBenders NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteBenders NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransBenders NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpBenders NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpBenders NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolBenders NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBenders)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, NULL, conshdlr, result, LP) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBenders)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, sol, conshdlr, result, RELAX) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBenders)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPconsBendersEnforceSolutions(scip, NULL, conshdlr, result, PSEUDO) );

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

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);

   (*result) = SCIP_FEASIBLE;
   performcheck = TRUE;
   infeasible = FALSE;


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
   if( performcheck )
   {
      for( i = 0; i < nactivebenders; i++ )
      {

         SCIP_CALL( SCIPsolveBendersSubproblems(scip, benders[i], sol, result, &infeasible, CHECK) );

         /* if the result is infeasible, it is not necessary to check any more subproblems. */
         if( (*result) == SCIP_INFEASIBLE )
            break;
      }

      /* in the case that the problem is feasible, this means that all subproblems are feasible. The auxiliary variables
       * still need to be updated. This is done by constructing a valid solution. */
      if( SCIPsolGetHeur(sol) == SCIPfindHeur(scip, "rins") )
         printf("result: %d, infeasible %d\n", *result, infeasible);
      if( (*result) == SCIP_FEASIBLE && infeasible )
      {
         if( !SCIPsolIsOriginal(sol) )
            SCIP_CALL( constructValidSolution(scip, conshdlr, sol) );
         (*result) = SCIP_INFEASIBLE;
      }
   }


   conshdlrdata->ncalls++;

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 0
static
SCIP_DECL_CONSPROP(consPropBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropBenders NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolBenders NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropBenders NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBenders)
{  /*lint --e{715}*/
   //SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   //SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveBenders NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveBenders NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableBenders NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableBenders NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsBenders NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintBenders NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyBenders NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseBenders NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsBenders NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsBenders NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsBenders)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of benders constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsBenders NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for benders constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBenders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             twophase            /**< should the two phase method be used? */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr;

   /* create benders constraint handler data */
   conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->ncalls = 0;

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
         conshdlrCopyBenders,
         consFreeBenders, consInitBenders, consExitBenders,
         consInitpreBenders, consExitpreBenders, consInitsolBenders, consExitsolBenders,
         consDeleteBenders, consTransBenders, consInitlpBenders,
         consSepalpBenders, consSepasolBenders, consEnfolpBenders, consEnforelaxBenders, consEnfopsBenders, consCheckBenders,
         consPropBenders, consPresolBenders, consRespropBenders, consLockBenders,
         consActiveBenders, consDeactiveBenders,
         consEnableBenders, consDisableBenders, consDelvarsBenders,
         consPrintBenders, consCopyBenders, consParseBenders,
         consGetVarsBenders, consGetNVarsBenders, consGetDiveBdChgsBenders, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBenders, consEnfopsBenders, consCheckBenders, consLockBenders,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveBenders) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBenders, consCopyBenders) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveBenders) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteBenders) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsBenders) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableBenders) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableBenders) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitBenders) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreBenders) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolBenders) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBenders) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsBenders) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsBenders) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsBenders) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitBenders) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreBenders) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolBenders) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpBenders) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseBenders) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolBenders, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintBenders) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropBenders, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropBenders) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpBenders, consSepasolBenders, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransBenders) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBenders) );

#endif

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
