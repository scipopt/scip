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
#define CONSHDLR_MAXSTALLS         1000 /**< disable after this many unsuccessful calls in a row */
#define DEFAULT_SOLBUFSIZE           10 /**< how many heuristic solutions should be buffered before we start checking */

static int ncurrentstalls = 0; /* number of times the exact lp was solved unsuccesfully in a row */

/** constraint handler data */
struct SolIntAssignment
{
   int*                  idx;                /**< variable indices that are fixed, sorted increasing */
   SCIP_Longint*         vals;               /**< values those vars were fixed to */
   int                   len;                /**< length of the two arrays */
};
typedef struct SolIntAssignment SOLINTASSIGNMENT;

struct SCIP_ConshdlrData
{
   SCIP_SOL**            solubuffer;         /**< buffer solutions for later checking here */
   SCIP_HASHTABLE*       solhash;            /**< hash solutions so we don't use the same integer assignment twice */
   SOLINTASSIGNMENT**    hashedassignments;  /**< array with all hashed assignments */
   int                   nhashedassignments; /**< number of elements in the array */
   int                   lenhash;            /**< length of the array */
   int                   lensolubuffer;      /**< size of solubuffer */
   int                   nbufferedsols;      /**< number of solutions currently in the solubuffer */
};

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeySol)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqSol)
{
   int i;
   SOLINTASSIGNMENT* sol1;
   SOLINTASSIGNMENT* sol2;

   assert(key1 != NULL);
   assert(key2 != NULL);

   sol1 = (SOLINTASSIGNMENT*)key1;
   sol2 = (SOLINTASSIGNMENT*)key2;

   if( sol1->len != sol2->len )
      return false;

   for( i = 0; i < sol1->len; i++ )
   {
      if( sol1->idx[i] != sol2->idx[i] || sol1->vals[i] != sol2->vals[i] )
         return FALSE;
   }

   return TRUE;
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(hashKeyValSol)
{
   SOLINTASSIGNMENT* sol;
   uint64_t signature;
   int i;

   sol = (SOLINTASSIGNMENT*)key;
   signature = 0;
   for( i = 0; i < sol->len; ++i )
      signature |= SCIPhashSignature64(sol->vals[i] * sol->idx[i]);

   return signature;
}

/** unlinks and copies a solution and adds it to the solution buffer */
static
SCIP_RETCODE bufferSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to add */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactsol constraint handler data */
   )
{
   SCIP_SOL* insertsol;

   SCIPdebugMessage("buffering solution from heuristic %s \n", SCIPheurGetName(SCIPsolGetHeur(sol)));

   SCIP_CALL( SCIPcreateSolCopy(scip, &insertsol, sol) );
   SCIP_CALL( SCIPunlinkSol(scip, insertsol) );

   /* extend solubuffer, if necessary */
   if( conshdlrdata->nbufferedsols == conshdlrdata->lensolubuffer )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->solubuffer, conshdlrdata->lensolubuffer, conshdlrdata->lensolubuffer * 2) );
      conshdlrdata->lensolubuffer *= 2;
   }

   /* put solution in buffer, according to objective value (best obj value goes last) */
   conshdlrdata->solubuffer[conshdlrdata->nbufferedsols] = insertsol;
   conshdlrdata->nbufferedsols++;

   return SCIP_OKAY;
}

/** frees all remaining solutions in buffer */
static
void clearSoluBuffer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactsol constraint handler data */
   )
{
   int i;

   for( i = 0; i < conshdlrdata->nbufferedsols; i++ )
   {
      SCIPfreeSol(scip, &conshdlrdata->solubuffer[i]);
   }

   conshdlrdata->nbufferedsols = 0;
}

/** creates assignment from integer variable-values in solution */
static
void solCreateSolAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to create assignment for */
   SOLINTASSIGNMENT**    assignment          /**< assignment array */
   )
{
   int nvars;
   int i;
   int solsize;
   SCIP_VAR** vars;

   assert(sol != NULL);
   assert(scip != NULL);

   SCIPallocBlockMemory(scip, assignment);
   SCIPallocClearBlockMemoryArray(scip, &(*assignment)->vals, SCIPgetNIntVars(scip) + SCIPgetNBinVars(scip));
   SCIPallocClearBlockMemoryArray(scip, &(*assignment)->idx, SCIPgetNIntVars(scip) + SCIPgetNBinVars(scip));

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   solsize = 0;

   for( i = 0; i < nvars; i++ )
   {
      if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_INTEGER && SCIPvarGetType(vars[i]) != SCIP_VARTYPE_BINARY )
         continue;

       (*assignment)->vals[solsize] = SCIPround(scip, SCIPgetSolVal(scip, sol, vars[i]));
       (*assignment)->idx[solsize] = SCIPvarGetIndex(vars[i]);
       solsize++;
   }

   (*assignment)->len = solsize;
}

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
   SCIP_CONS** consprob;
   SCIP_SOL* exactsol;
   SOLINTASSIGNMENT* assignment;
   SCIP_SOL* worksol;
   SCIP_Bool foundsol;
   SCIP_Bool lperror;
   int nintvars;
   int nvars;
   int nfixedvars;
   int nconsprob;
   int i;
   int c;
   SCIP_Real starttime;
   SCIP_Real endtime;
   SCIP_CONSHDLRDATA* conshdlrdata;
#ifdef NDEBUG
   SCIP_RETCODE retstat;
#endif

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;
   foundsol = FALSE;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   assert(conshdlrdata != NULL);

   /* disable exact sol if we stalled too often in a row */
   if( ncurrentstalls >= CONSHDLR_MAXSTALLS )
      return SCIP_OKAY;

   /* if the solution doesn't come from a heuristic, ignore it */
   if( SCIPsolGetType(sol) != SCIP_SOLTYPE_HEUR )
      return SCIP_OKAY;

   if( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   if( SCIPgetNContVars(scip) > 0.8 * SCIPgetNVars(scip) )
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

   /** if we are at a point where we can't dive, buffer the solution and return */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || !SCIPtreeHasCurrentNodeLP(scip->tree) || SCIPnodeGetType(SCIPgetCurrentNode(scip)) != SCIP_NODETYPE_FOCUSNODE )
   {
      SCIP_CALL( bufferSolution(scip, sol, conshdlrdata) );
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   if( !SCIPisLPConstructed(scip) )
   {
      SCIP_Bool cutoff = FALSE;
      SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
      SCIP_CALL( SCIPflushLP(scip) );
   }

   /** no point trying to repair solutions that are already exact */
   if( SCIPisExactSol(scip, sol) )
      return SCIP_OKAY;

   starttime = SCIPgetSolvingTime(scip);
   nconsprob = SCIPgetNConss(scip);
   consprob = SCIPgetConss(scip);

   /* check if solution is floating point feasible */
   for( c = 0; c < nconsprob && *result == SCIP_FEASIBLE; ++c )
   {
      SCIP_Real activity;
      SCIP_ROW* row;

      /* get row corresponding to constraint */
      row = SCIPconsGetRow(scip, consprob[c]);
      assert(row != NULL);

      /* get row activity */
      activity = SCIPgetRowSolActivity(scip, row, sol);

      /* check if the constraint is violated */
      if( SCIPisFeasLT(scip, activity, SCIProwGetLhs(row)) || SCIPisFeasGT(scip, activity, SCIProwGetRhs(row)) )
         *result = SCIP_INFEASIBLE;
   }

   /* do not check infeasible solutions */
   if( *result == SCIP_INFEASIBLE )
      return SCIP_OKAY;

   /* first, check if we already tried a solution with this integer assignment */
   solCreateSolAssignment(scip, sol, &assignment);
   if( SCIPhashtableExists(conshdlrdata->solhash, (void*) assignment) )
   {
      SCIPdebugMessage("rejecting solution that was already checked \n");
      SCIPdebug(SCIPprintSol(scip, sol, NULL, 0));

      SCIPfreeBlockMemoryArray(scip, &assignment->idx, SCIPgetNIntVars(scip));
      SCIPfreeBlockMemoryArray(scip, &assignment->vals, SCIPgetNIntVars(scip));
      SCIPfreeBlockMemory(scip, &assignment);
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }
   else
   {
      SCIPdebugMessage("checking solution for the first time: \n");
      SCIPdebug(SCIPprintSol(scip, sol, NULL, 0));
      if( conshdlrdata->lenhash == conshdlrdata->nhashedassignments )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->hashedassignments, conshdlrdata->lenhash, conshdlrdata->lenhash * 2) );
         conshdlrdata->lenhash *= 2;
      }
      conshdlrdata->hashedassignments[conshdlrdata->nhashedassignments] = assignment;
      conshdlrdata->nhashedassignments++;

      SCIPhashtableInsert(conshdlrdata->solhash, assignment);
   }

   /* next, check if we should buffer the solution instead of checking it right now */
   if( conshdlrdata->nbufferedsols < DEFAULT_SOLBUFSIZE )
   {
      SCIP_CALL( bufferSolution(scip, sol, conshdlrdata) );
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   scip->stat->ncallsexactsol++;

   /* start exact diving */
   SCIP_CALL( SCIPstartExactDive(scip) );

   /* sort solubuffer by objva; try to repair best solutions first */
   SCIPsortPtr((void**)conshdlrdata->solubuffer, SCIPsolComp, conshdlrdata->nbufferedsols);

   while( conshdlrdata->nbufferedsols > 0 && !foundsol )
   {
      /* best solution is last in solubuffer */
      worksol = conshdlrdata->solubuffer[conshdlrdata->nbufferedsols - 1];

      /* only try to fix solutions that improve the cutoffbound */
      if( !SCIPisLT(scip, SCIPgetSolTransObj(scip, worksol), SCIPgetUpperbound(scip)) )
      {
         SCIPdebugMessage("don't repair heuristic with obj value (%g) greater than upper bound (%g) \n",  SCIPgetSolTransObj(scip, worksol), SCIPgetUpperbound(scip));
         SCIP_CALL( SCIPfreeSol(scip, &worksol) );
         conshdlrdata->nbufferedsols--;
         continue;
      }

      SCIPdebugMessage("attempting to repair solution from heuristic %s with floating point objval %g \n", SCIPheurGetName(SCIPsolGetHeur(worksol)), SCIPgetSolTransObj(scip, worksol));

      /* set the bounds of the variables: fixed for integers, global bounds for continuous */
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);
      nintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
      nfixedvars = 0;

      for( i = 0; i < nintvars; ++i )
      {
         if( SCIPvarGetStatusExact(vars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_Real solval;
            solval = SCIPgetSolVal(scip, worksol, vars[i]);
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

               nfixedvars++;

               RatFreeBuffer(SCIPbuffer(scip), &newbound);
            }
            else
               *result = SCIP_INFEASIBLE;
         }
      }

      if( *result == SCIP_INFEASIBLE )
      {
         SCIP_CALL( SCIPfreeSol(scip, &worksol) );
         conshdlrdata->nbufferedsols--;
         continue;
      }

      *result = SCIP_INFEASIBLE;

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

         SCIPsolSetHeur(exactsol, SCIPsolGetHeur(worksol));
         SCIP_CALL( SCIPtrySolFreeExact(scip, &exactsol, FALSE, FALSE, FALSE, FALSE, TRUE, &foundsol) );

         /* if we found a solution we do not try to complete the others, since they have worse objective values */
         if( foundsol )
         {
            clearSoluBuffer(scip, conshdlrdata);
         }
         else
         {
            SCIPfreeSol(scip, &worksol);
            conshdlrdata->nbufferedsols--;
         }
      }
   }

   /* terminate exact diving */
   SCIP_CALL( SCIPendExactDive(scip) );

   endtime = SCIPgetSolvingTime(scip);
   if( foundsol )
   {
      SCIPdebugMessage("successfully found feasible improving solution, objval %g, upperbound %g\n", SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip));
      ncurrentstalls = 0;
      scip->stat->nfoundexactsol++;
      scip->stat->timesuccessexactsol += (endtime - starttime);
   }
   else
   {
      SCIPdebugMessage("repaired solution not feasible or not improving, objval %g, upperbound %g \n", SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip));
      ncurrentstalls++;
      scip->stat->timefailexactsol += (endtime - starttime);
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExactSol)
{  /*lint --e{715}*/

   /* do nothing since we are not handling constraints */
   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
/**! [SnippetConsFreeKnapsack] */
static
SCIP_DECL_CONSFREE(consFreeExactSol)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitExactSol)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nvars;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create hashdata for integer assignments */
   SCIPallocBlockMemoryArray(scip, &conshdlrdata->hashedassignments, 10);
   SCIP_CALL( SCIPhashtableCreate(&(conshdlrdata->solhash), SCIPblkmem(scip), 10, hashGetKeySol, hashKeyEqSol, hashKeyValSol, NULL) );

   conshdlrdata->nhashedassignments = 0;
   conshdlrdata->lenhash = 10;

   /* allocate data for solution buffer */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->solubuffer, DEFAULT_SOLBUFSIZE) );
   conshdlrdata->lensolubuffer = DEFAULT_SOLBUFSIZE;
   conshdlrdata->nbufferedsols = 0;

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitExactSol)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free solution hashdata */
   SCIPhashtableRemoveAll(conshdlrdata->solhash);
   SCIPhashtableFree(&(conshdlrdata->solhash));
   for( i = 0; i < conshdlrdata->nhashedassignments; i++ )
   {
      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->hashedassignments[i]->idx, conshdlrdata->hashedassignments[i]->len);
      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->hashedassignments[i]->vals, conshdlrdata->hashedassignments[i]->len);
      SCIPfreeBlockMemory(scip, &conshdlrdata->hashedassignments[i]);
   }
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->hashedassignments, conshdlrdata->lenhash);
   conshdlrdata->nhashedassignments = 0;

   /* free solubuffer */
   for( i = 0; i < conshdlrdata->nbufferedsols; i++ )
   {
      SCIPfreeSol(scip, &conshdlrdata->solubuffer[i]);
   }
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->solubuffer, conshdlrdata->lensolubuffer);
   conshdlrdata->nbufferedsols = 0;

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
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpExactSol, consEnfopsExactSol, consCheckExactSol, consLockExactSol,
         conshdlrdata) );

   SCIPconshdlrSetInit(conshdlr, consInitExactSol);
   SCIPconshdlrSetExit(conshdlr, consExitExactSol);
   SCIPconshdlrSetFree(conshdlr, consFreeExactSol);
   assert(conshdlr != NULL);

   return SCIP_OKAY;
}
