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
#include "scip/cons.h"
#include "scip/scip_exact.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_sol.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/stat.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "exactsol"
#define CONSHDLR_DESC          "constraint to ensure that primal solutions report back exact solutions"
#define CONSHDLR_ENFOPRIORITY  -9999999 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */


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
   SCIP_SOL* exactsol;
   SCIP_Bool foundsol;
   SCIP_Bool lperror;
   int nintvars;
   int nvars;
   int i;
   int c;
   SCIP_Real starttime;
   SCIP_Real endtime;
#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;
   foundsol = FALSE;

   /* if the solution doesn't come from a heuristic, ignore it */
   if( SCIPsolGetType(sol) != SCIP_SOLTYPE_HEUR )
      return SCIP_OKAY;

   /* if we are not solving exactly, we have nothing to check */
   if( !SCIPisExactSolve(scip) )
      return SCIP_OKAY;

   /* if we're already in exact diving mode, we already computed an exact solution with this constraint handler and
    * are checking if it's actually feasible */
   if( SCIPlpExactDiving(scip->lpexact) )
      return SCIP_OKAY;

   /* we also don't want to execute the handler, if we are in "normal" diving mode */
   if( SCIPlpDiving(scip->lp) )
      return SCIP_OKAY;

   /** @todo exip: buffer solutions in this case instead */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_Bool cutoff = FALSE;

      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
      SCIP_CALL( SCIPflushLP(scip) );
   }

   if( !SCIPtreeHasCurrentNodeLP(scip->tree))
      return SCIP_OKAY;

   if( !SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_FOCUSNODE )
      return SCIP_OKAY;

   if( SCIPisExactSol(scip, sol) )
      return SCIP_OKAY;

   starttime = SCIPgetSolvingTime(scip);

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
   scip->stat->ncallsexactsol++;

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

            SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], SCIPround(scip, solval)) );
            SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], SCIPround(scip, solval)) );

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

   /** @todo exip handle the ray case to be included */

   /* check if this is a feasible solution */
   if( !lperror && SCIPgetLPExactSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
   {
      /** @todo exip: replace these 2 methods with proper wrapped ones from scip_sol.h */
      SCIP_CALL( SCIPsolCreateLPSolExact(&exactsol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
          scip->tree, scip->lpexact, NULL) );
      SCIP_CALL( SCIPsolOverwriteFPSolWithExact(exactsol, scip->set, scip->stat, scip->origprob, scip->transprob, scip->tree) );

      SCIPsolSetHeur(exactsol, SCIPsolGetHeur(sol));
      SCIP_CALL( SCIPtrySolFreeExact(scip, &exactsol, FALSE, FALSE, FALSE, FALSE, TRUE, &foundsol) );

      if( foundsol )
      {
         *result = SCIP_FEASIBLE;
      }
   }

   /* terminate exact diving */
   SCIP_CALL( SCIPendExactDive(scip) );

   endtime = SCIPgetSolvingTime(scip);
   if( foundsol )
   {
      scip->stat->nfoundexactsol++;
      scip->stat->timesuccessexactsol += (endtime - starttime);
   }
   else
      scip->stat->timefailexactsol += (endtime - starttime);

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
   conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpExactSol, consEnfopsExactSol, consCheckExactSol, consLockExactSol,
         conshdlrdata) );
   assert(conshdlr != NULL);

   return SCIP_OKAY;
}