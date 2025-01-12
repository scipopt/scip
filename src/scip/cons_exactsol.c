/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_exactsol.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for ensuring that primal solution is exact
 * @author Antonia Chmiela
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_exactsol.h"
#include "scip/cons_exactlp.h"
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
#define DEFAULT_CHECKFPFEASIBILITY TRUE /**< should a solution be checked in floating-point arithmetic prior to being processed? */

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
   int                   nhashedassignments; /**< number of elements in the hashedassignments array */
   int                   lenhash;            /**< length of the hashedassignments array */
   int                   nbufferedsols;      /**< number of solutions currently in the solubuffer */
   int                   lensolubuffer;      /**< length of the solubuffer */
   SCIP_Bool             checkfpfeasibility; /**< should a solution be checked in floating-point arithmetic prior to being processed? */
   int                   probhasequations;   /**< does the problem have equations? (-1 unknown, 0 no, 1 yes) */
};

/** gets the key of the given element */
static
SCIP_DECL_HASHGETKEY(hashGetKeyAssignment)
{  /*lint --e{715}*/
   /* the key is the element itself */
   return elem;
}

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(hashKeyEqAssignment)
{ /*lint --e{715}*/
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
SCIP_DECL_HASHKEYVAL(hashKeyValAssignment)
{ /*lint --e{715}*/
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
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->solubuffer,
            conshdlrdata->lensolubuffer, conshdlrdata->lensolubuffer * 2) );
      conshdlrdata->lensolubuffer *= 2;
   }

   /* put solution in buffer */
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
      SCIP_CALL_ABORT( SCIPfreeSol(scip, &conshdlrdata->solubuffer[i]) );
   }

   conshdlrdata->nbufferedsols = 0;
}

/** creates assignment from integer variable-values in solution */
static
void solCreateSolAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to create assignment for */
   SOLINTASSIGNMENT**    assignment          /**< address of assignment */
   )
{ /*lint --e{522, 776}*/
   int nvars;
   int i;
   int solsize;
   SCIP_VAR** vars;

   assert(sol != NULL);
   assert(scip != NULL);

   SCIPallocBlockMemory(scip, assignment);
   if( *assignment == NULL )
      return;
   SCIPallocClearBlockMemoryArray(scip, &(*assignment)->vals, SCIPgetNIntVars(scip) + SCIPgetNBinVars(scip));
   SCIPallocClearBlockMemoryArray(scip, &(*assignment)->idx, SCIPgetNIntVars(scip) + SCIPgetNBinVars(scip));

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   solsize = 0;

   for( i = 0; i < nvars; i++ )
   {
      if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_INTEGER && SCIPvarGetType(vars[i]) != SCIP_VARTYPE_BINARY )
         continue;

      (*assignment)->vals[solsize] = (SCIP_Longint) SCIPround(scip, SCIPgetSolVal(scip, sol, vars[i]));
      (*assignment)->idx[solsize] = SCIPvarGetIndex(vars[i]);
      solsize++;
   }

   (*assignment)->len = solsize;
   assert(solsize == SCIPgetNIntVars(scip) + SCIPgetNBinVars(scip));
}

/** creates assignment from integer variable-values in solution */
static
void solFreeAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   SOLINTASSIGNMENT**    assignment          /**< address of assignment */
   )
{
   assert(scip != NULL);
   assert(*assignment != NULL);

   SCIPfreeBlockMemoryArray(scip, &(*assignment)->idx, (*assignment)->len);
   SCIPfreeBlockMemoryArray(scip, &(*assignment)->vals, (*assignment)->len);
   SCIPfreeBlockMemory(scip, assignment);
}

static
void setProbHasEquations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< exactsol constraint handler data */
   )
{
   SCIP_CONS** conss;
   int nconss;
   int c;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->probhasequations != -1 )
      return;

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   success = TRUE;

   conshdlrdata->probhasequations = 0;

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsGetHdlr(conss[c]) == SCIPfindConshdlr(scip, "linear-exact") )
      {
         // constraint is an equality constraint
         if( RatIsEqual(SCIPconsGetRhsExact(scip, conss[c], &success), SCIPconsGetLhsExact(scip, conss[c], &success)) ) /*lint !e864*/
         {
            // check if there are non-integer variables
            SCIP_VAR** vars = SCIPgetVarsExactLinear(scip, conss[c]);

            for( int i = 0; i < SCIPgetNVarsExactLinear(scip, conss[c]); ++i )
            {
               if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_INTEGER && SCIPvarGetType(vars[i]) != SCIP_VARTYPE_BINARY )
               {
                  conshdlrdata->probhasequations = 1;
                  break;
               }
            }
         }
         conshdlrdata->probhasequations = 1;
         break;
      }
      else
      {
         // unspported constraint type -> throw error
         SCIPerrorMessage("Unsupported constraint type in exactsol constraint handler: %s\n", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])));
         SCIPABORT();
      }
   }
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

/** constraint enforcing method of constraint handler for LP solutions */
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
   SOLINTASSIGNMENT* assignment = NULL;
   SCIP_SOL* worksol;
   SCIP_Bool foundsol;
   SCIP_Bool lperror;
   SCIP_Bool checkfpfeasibility;
   int nintvars;
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

   /* if we are not solving exactly, we have nothing to check */
   if( !SCIPisExactSolve(scip) )
      return SCIP_OKAY;

   if( conshdlrdata->probhasequations == -1 )
      setProbHasEquations(scip, conshdlrdata);

   if( conshdlrdata->probhasequations == 1 )
      return SCIP_OKAY;

   /* disable exact sol if we stalled too often in a row */
   if( ncurrentstalls >= CONSHDLR_MAXSTALLS )
      return SCIP_OKAY;

   /* if the solution doesn't come from a heuristic, ignore it */
   if( SCIPsolGetType(sol) != SCIP_SOLTYPE_HEUR )
      return SCIP_OKAY;

   /* do not run for problems that contain mostly continuous variables */
   if( SCIPgetNContVars(scip) > 0.8 * SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* do not run for problems that are purely integer */
   if( SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

   /* if we're already in exact diving mode, we already computed an exact solution
    * with this constraint handler and are checking if it's actually feasible */
   if( SCIPlpExactDiving(scip->lpexact) )
      return SCIP_OKAY;

   /* we also don't want to execute the handler, if we are in "normal" diving mode */
   if( SCIPlpDiving(scip->lp) )
      return SCIP_OKAY;

   /* do not run for trivial since it would only be buffered and never used */
   if( strcmp(SCIPheurGetName(SCIPsolGetHeur(sol)), "trivial") == 0 )
      return SCIP_OKAY;

   /* do not run for solutions that are already exact */
   if( SCIPisExactSol(scip, sol) )
      return SCIP_OKAY;

   /* if we are at a point where we can't dive exactly, buffer the solution and return */
   if( !SCIPtreeHasCurrentNodeLP(scip->tree) || SCIPnodeGetType(SCIPgetCurrentNode(scip)) != SCIP_NODETYPE_FOCUSNODE )
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

   /* no point trying to repair solutions that are already exact */
   if( SCIPisExactSol(scip, sol) )
      return SCIP_OKAY;

   starttime = SCIPgetSolvingTime(scip);
   nconsprob = SCIPgetNConss(scip);
   consprob = SCIPgetConss(scip);

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/" CONSHDLR_NAME "/checkfpfeasibility", &checkfpfeasibility) );

   /* check if solution is floating point feasible */
   for( c = 0; c < nconsprob && *result == SCIP_FEASIBLE && checkfpfeasibility; ++c )
   {
      SCIP_Real activity;
      SCIP_ROW* row;

      /* get row corresponding to constraint */
      row = SCIPconsGetRow(scip, consprob[c]);
      if( row == NULL )
         continue;

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
   if( assignment != NULL && SCIPhashtableExists(conshdlrdata->solhash, (void*) assignment) )
   {
      SCIPdebugMessage("rejecting solution that was already checked \n");
      SCIPdebug(SCIPprintSol(scip, sol, NULL, 0));

      solFreeAssignment(scip, &assignment);
      *result = SCIP_INFEASIBLE;

      return SCIP_OKAY;
   }
   else
   {
      SCIPdebugMessage("checking solution for the first time: \n");
      SCIPdebug(SCIPprintSol(scip, sol, NULL, 0));
      /* add assignment to the hashtable, extend assignment array, if necessary */
      if( conshdlrdata->lenhash == conshdlrdata->nhashedassignments )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &conshdlrdata->hashedassignments,
               conshdlrdata->lenhash, conshdlrdata->lenhash * 2) );
         conshdlrdata->lenhash *= 2;
      }
      conshdlrdata->hashedassignments[conshdlrdata->nhashedassignments] = assignment;
      conshdlrdata->nhashedassignments++;

      SCIP_CALL( SCIPhashtableInsert(conshdlrdata->solhash, assignment) );
   }

   /* next, check if we should buffer the solution instead of checking it right now, done in one of the following cases
      - buffer is not full and solution not more than 20% improving
      - exact diving not possible at this point in time (mostly if lp state is not clean) */
   if( conshdlrdata->nbufferedsols < DEFAULT_SOLBUFSIZE || !SCIPisExactDivePossible(scip) )
   {
      SCIP_Real multiplier;

      SCIP_CALL( bufferSolution(scip, sol, conshdlrdata) );

      multiplier = SCIPgetSolTransObj(scip, sol) > 0 ? 1.2 : 0.8;
      /* if the new solution is at least 20% better than the current upperbound, we stop buffering and repair immediately */
      if( !SCIPisLT(scip, multiplier * SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip)) || !SCIPisExactDivePossible(scip) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }

   scip->stat->ncallsexactsol++;

   /* start exact diving */
   SCIP_CALL( SCIPstartExactDive(scip) );

   /* sort solubuffer by objval try to repair best solutions first */
   SCIPsortPtr((void**)conshdlrdata->solubuffer, SCIPsolComp, conshdlrdata->nbufferedsols);

   while( conshdlrdata->nbufferedsols > 0 && !foundsol )
   {
      /* best solution is last in solubuffer */
      worksol = conshdlrdata->solubuffer[conshdlrdata->nbufferedsols - 1];

      /* only try to fix solutions that improve the cutoffbound */
      if( !SCIPisLT(scip, SCIPgetSolTransObj(scip, worksol), SCIPgetUpperbound(scip)) )
      {
         SCIPdebugMessage("don't repair heuristic with obj value (%g) greater than upper bound (%g) \n",
            SCIPgetSolTransObj(scip, worksol), SCIPgetUpperbound(scip));
         SCIP_CALL( SCIPfreeSol(scip, &worksol) );
         conshdlrdata->nbufferedsols--;
         continue;
      }

      SCIPdebugMessage("attempting to repair solution from heuristic %s with floating point objval %g \n",
         SCIPheurGetName(SCIPsolGetHeur(worksol)), SCIPgetSolTransObj(scip, worksol));

      /* set the bounds of the variables: fixed for integers, global bounds for continuous */
      vars = SCIPgetVars(scip);
      nintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

      for( i = 0; i < nintvars; ++i )
      {
         if( SCIPvarGetStatusExact(vars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_Real solval;
            solval = SCIPgetSolVal(scip, worksol, vars[i]);

            assert(RatIsLE(SCIPvarGetLbLocalExact(vars[i]), SCIPvarGetUbLocalExact(vars[i])));

            /* check if solution value is supposed to be integral */
            if( SCIPisIntegral(scip, solval) )
            {
               SCIP_Rational* newbound;

               SCIP_CALL( RatCreateBuffer(SCIPbuffer(scip), &newbound) );

               /* create rational solval and round it to the nearest integer */
               RatSetReal(newbound, solval);
               RatRound(newbound, newbound, SCIP_R_ROUND_NEAREST);

               SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], SCIPround(scip, solval)) );
               SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], SCIPround(scip, solval)) );

               SCIP_CALL( SCIPchgVarLbExactDive(scip, vars[i], newbound) );
               SCIP_CALL( SCIPchgVarUbExactDive(scip, vars[i], newbound) );

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
         SCIP_CALL( SCIPcreateLPSolExact(scip, &exactsol, NULL) );
         SCIP_CALL( SCIPoverwriteFPsol(scip, exactsol) );

         SCIPsolSetHeur(exactsol, SCIPsolGetHeur(worksol));
         SCIP_CALL( SCIPtrySolFreeExact(scip, &exactsol, FALSE, FALSE, FALSE, FALSE, TRUE, &foundsol) );

         /* if we found a solution we do not try to complete the others, since they have worse objective values */
         if( foundsol )
         {
            clearSoluBuffer(scip, conshdlrdata);
         }
         else
         {
            SCIP_CALL( SCIPfreeSol(scip, &worksol) );
            conshdlrdata->nbufferedsols--;
         }
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &worksol) );
         conshdlrdata->nbufferedsols--;
      }
   }

   /* terminate exact diving */
   SCIP_CALL( SCIPendExactDive(scip) );

   endtime = SCIPgetSolvingTime(scip);
   if( foundsol )
   {
      SCIPdebugMessage("successfully found feasible improving solution, objval %g, upperbound %g\n",
         SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip));
      ncurrentstalls = 0;
      scip->stat->nfoundexactsol++;
      scip->stat->timesuccessexactsol += (endtime - starttime);
   }
   else
   {
      SCIPdebugMessage("repaired solution not feasible or not improving, objval %g, upperbound %g \n",
         SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip));
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
{  /*lint --e{715, 522}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* create hashdata for integer assignments */
   SCIPallocBlockMemoryArray(scip, &conshdlrdata->hashedassignments, DEFAULT_SOLBUFSIZE);
   SCIP_CALL( SCIPhashtableCreate(&(conshdlrdata->solhash), SCIPblkmem(scip), DEFAULT_SOLBUFSIZE, hashGetKeyAssignment, hashKeyEqAssignment, hashKeyValAssignment, NULL) );

   conshdlrdata->nhashedassignments = 0;
   conshdlrdata->lenhash = DEFAULT_SOLBUFSIZE;

   /* allocate data for solution buffer */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->solubuffer, DEFAULT_SOLBUFSIZE) );
   conshdlrdata->lensolubuffer = DEFAULT_SOLBUFSIZE;
   conshdlrdata->nbufferedsols = 0;
   conshdlrdata->probhasequations = -1;

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitExactSol)
{  /*lint --e{715, 866}*/
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
      SCIP_CALL( SCIPfreeSol(scip, &conshdlrdata->solubuffer[i]) );
   }
   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->solubuffer, conshdlrdata->lensolubuffer);
   conshdlrdata->nbufferedsols = 0;

   return SCIP_OKAY;
}

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyExactSol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrExactSol(scip) );

   *valid = TRUE;

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

   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyExactSol, NULL) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxExactSol) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/checkfpfeasibility",
         "should a solution be checked in floating-point arithmetic prior to being processed?",
         &conshdlrdata->checkfpfeasibility, TRUE, DEFAULT_CHECKFPFEASIBILITY, NULL, NULL) );

   SCIPconshdlrSetInit(conshdlr, consInitExactSol);
   SCIPconshdlrSetExit(conshdlr, consExitExactSol);
   SCIPconshdlrSetFree(conshdlr, consFreeExactSol);
   assert(conshdlr != NULL);

   return SCIP_OKAY;
}
