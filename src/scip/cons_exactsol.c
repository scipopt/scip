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

#include "scip/def.h"
#include "scip/cons_exactlinear.h"
#include "scip/cons_exactsol.h"
#include "scip/pub_cons.h"
#include "scip/pub_heur.h"
#include "scip/pub_lp.h"
#include "scip/pub_lpexact.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/rational.h"
#include "scip/scip_cons.h"
#include "scip/scip_exact.h"
#include "scip/scip_lp.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/set.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "exactsol"
#define CONSHDLR_DESC          "constraint handler for repairing floating-point primal solutions to satisfy exact feasibility"
#define CONSHDLR_ENFOPRIORITY  -9999999 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_CHECKFPFEASIBILITY TRUE /**< should a solution be checked in floating-point arithmetic prior to being processed? */
#define DEFAULT_MAXSTALLS          1000 /**< maximal number of consecutive repair calls without success */
#define DEFAULT_SOLBUFSIZE           10 /**< size of solution buffer */
#define DEFAULT_MINIMPROVE          0.2 /**< minimal percentage of primal improvement to trigger solution processing */

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
   int                   probhasconteqs;     /**< does the problem have equations with continuous variables? (-1 unknown, 0 no, 1 yes) */
   int                   ncurrentstalls;     /**< number of times the exact lp was solved unsuccessfully in a row */
   SCIP_Bool             checkfpfeasibility; /**< should a solution be checked in floating-point arithmetic prior to being processed? */
   int                   maxstalls;          /**< maximal number of consecutive repair calls without success */
   int                   solbufsize;         /**< size of solution buffer */
   SCIP_Real             minimprove;         /**< minimal percentage of primal improvement to trigger solution processing */
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
   SOLINTASSIGNMENT* sol1;
   SOLINTASSIGNMENT* sol2;
   int i;

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
SCIP_RETCODE solCreateSolAssignment(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to create assignment for */
   SOLINTASSIGNMENT**    assignment          /**< address of assignment */
   )
{ /*lint --e{522, 776}*/
   SCIP_VAR** vars;
   int nvars;
   int ncontvars;
   int nintegers;
   int i;

   assert(sol != NULL);
   assert(scip != NULL);

   /* get all problem variables and integer region in vars array */
   SCIP_CALL( SCIPgetSolVarsData(scip, sol, &vars, &nvars, NULL, NULL, NULL, NULL, NULL, &ncontvars) );
   nintegers = nvars - ncontvars;
   assert(nintegers >= 0);

   SCIP_CALL( SCIPallocBlockMemory(scip, assignment) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(*assignment)->vals, nintegers) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(*assignment)->idx, nintegers) );

   for( i = 0; i < nintegers; i++ )
   {
      assert(SCIPvarIsIntegral(vars[i]));

      (*assignment)->vals[i] = (SCIP_Longint) SCIPround(scip, SCIPgetSolVal(scip, sol, vars[i]));
      (*assignment)->idx[i] = SCIPvarGetIndex(vars[i]);
   }

   (*assignment)->len = nintegers;

   return SCIP_OKAY;
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

/** checks whether equation constraints with non-integral variables are present */
static
void checkProbHasContEqs(
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

   if( conshdlrdata->probhasconteqs != -1 )
      return;

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   success = TRUE;

   conshdlrdata->probhasconteqs = 0;

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsGetHdlr(conss[c]) == SCIPfindConshdlr(scip, "exactlinear") )
      {
         /* constraint is an equality constraint */
         if( SCIPrationalIsEQ(SCIPconsGetRhsExact(scip, conss[c], &success), SCIPconsGetLhsExact(scip, conss[c], &success)) ) /*lint !e864*/
         {
            /* check if there are continuous variables involved */
            SCIP_VAR** vars = SCIPgetVarsExactLinear(scip, conss[c]);
            int nvars = SCIPgetNVarsExactLinear(scip, conss[c]);

            for( int i = 0; i < nvars; ++i )
            {
               if( !SCIPvarIsIntegral(vars[i]) )
               {
                  conshdlrdata->probhasconteqs = 1;
                  break;
               }
            }
         }
         conshdlrdata->probhasconteqs = 1;
         break;
      }
      else
      {
         /* unsupported constraint type -> throw error */
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
   int nvars;
   int nintegers;
   int nconsprob;
   int i;
   int c;
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
   if( !SCIPisExact(scip) )
      return SCIP_OKAY;

   /**@todo add event handler to check again if constraints were added/modified or a variable (impl) type changed */
   if( conshdlrdata->probhasconteqs == -1 )
      checkProbHasContEqs(scip, conshdlrdata);

   /* disable exact sol if we stalled too often in a row */
   if( conshdlrdata->ncurrentstalls >= conshdlrdata->maxstalls )
      return SCIP_OKAY;

   /* if the solution doesn't come from a heuristic, ignore it */
   if( SCIPsolGetType(sol) != SCIP_SOLTYPE_HEUR )
      return SCIP_OKAY;

   /* do not run if the solution comes from the trivial heuristic for the following reason: it typically creates the
    * first solution, which would be processed immediately, because it improves the primal bound by an infinite amount;
    * however, its quality is usually bad and superseeded quickly by solutions from other heuristics
    */
   if( strcmp(SCIPheurGetName(SCIPsolGetHeur(sol)), "trivial") == 0 )
      return SCIP_OKAY;

   /* do not run for problems that contain mostly continuous variables */
   if( SCIPgetNContVars(scip) > 0.8 * SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* do not run for problems that are purely integer */
   if( SCIPgetNContVars(scip) == 0 )
      return SCIP_OKAY;

   /* if we're already in exact diving mode, we already computed an exact solution
    * with this constraint handler and are checking if it's actually feasible
    */
   if( SCIPinExactDive(scip) )
      return SCIP_OKAY;

   /* we also don't want to execute the handler, if we are in "normal" diving mode */
   if( SCIPinDive(scip) )
      return SCIP_OKAY;

   /* do not run for solutions that are already exact */
   if( SCIPisExactSol(scip, sol) )
      return SCIP_OKAY;

   /* if we are at a point where we can't dive exactly, buffer the solution and return */
   if( !SCIPhasCurrentNodeLP(scip) || SCIPnodeGetType(SCIPgetCurrentNode(scip)) != SCIP_NODETYPE_FOCUSNODE )
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

   nconsprob = SCIPgetNConss(scip);
   consprob = SCIPgetConss(scip);

   /* check if solution is floating-point feasible */
   if( conshdlrdata->checkfpfeasibility )
   {
      for( c = 0; c < nconsprob && *result == SCIP_FEASIBLE ; ++c )
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

      /* do not continue for floating-point infeasible solutions */
      if( *result == SCIP_INFEASIBLE )
         return SCIP_OKAY;
   }

   /* first, check if we already tried a solution with this integer assignment */
   SCIP_CALL( solCreateSolAssignment(scip, sol, &assignment) );
   if( assignment != NULL && SCIPhashtableExists(conshdlrdata->solhash, (void*) assignment) )
   {
      SCIPdebugMessage("rejecting solution that was already checked\n");
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

   /* add solution to buffer */
   SCIP_CALL( bufferSolution(scip, sol, conshdlrdata) );

   /* stop if exact diving is not possible at this point in time (mostly if lp state is not clean) */
   if( !SCIPisExactDivePossible(scip) )
   {
      *result = SCIP_INFEASIBLE;
      return SCIP_OKAY;
   }

   /* stop if the new solution does not improve the current upperbound sufficiently and the buffer is not full;
    * otherwise we continue by processing the buffer
    */
   if( conshdlrdata->nbufferedsols < DEFAULT_SOLBUFSIZE )
   {
      SCIP_Real multiplier;

      multiplier = SCIPgetSolTransObj(scip, sol) > 0 ? 1 + conshdlrdata->minimprove : 1 - conshdlrdata->minimprove;
      if( !SCIPisLT(scip, multiplier * SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip)) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }

   /* start exact diving and set global bounds of continuous variables */
   SCIP_CALL( SCIPstartExactDive(scip) );

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nintegers = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   for( i = nintegers; i < nvars; ++i )
   {
      if( SCIPvarGetStatusExact(vars[i]) == SCIP_VARSTATUS_COLUMN )
      {
         SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], SCIPvarGetLbGlobal(vars[i])) );
         SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], SCIPvarGetUbGlobal(vars[i])) );

         SCIP_CALL( SCIPchgVarLbExactDive(scip, vars[i], SCIPvarGetLbGlobalExact(vars[i])) );
         SCIP_CALL( SCIPchgVarUbExactDive(scip, vars[i], SCIPvarGetUbGlobalExact(vars[i])) );
      }
   }

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

      /* set the bounds of the variables: fixed for integral variables, global bounds for continuous ones */
      for( i = 0; i < nintegers; ++i )
      {
         if( SCIPvarGetStatusExact(vars[i]) == SCIP_VARSTATUS_COLUMN )
         {
            SCIP_Real solval;
            solval = SCIPgetSolVal(scip, worksol, vars[i]);

            assert(SCIPrationalIsLE(SCIPvarGetLbLocalExact(vars[i]), SCIPvarGetUbLocalExact(vars[i])));

            /* check if solution value is integral and abort if not, except if integrality is weakly implied: then the
             * solution value could be fractional in a floating-point feasible solution and we know that an optimal
             * solution with integral value exist; in this case we currently round and fix its value
             */
            /**@todo once implied integrality detection is made exact, test whether it improves performance to leave
             *       continuous implied integral variables unfixed or fix them only if they take a nearly integral value
             */
            if( SCIPisIntegral(scip, solval) || SCIPvarGetImplType(vars[i]) == SCIP_IMPLINTTYPE_WEAK )
            {
               SCIP_RATIONAL* newbound;

               SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &newbound) );

               /* create rational solval and round it to the nearest integer */
               SCIPrationalSetReal(newbound, solval);
               SCIPrationalRoundInteger(newbound, newbound, SCIP_R_ROUND_NEAREST);

               SCIP_CALL( SCIPchgVarLbDive(scip, vars[i], SCIPround(scip, solval)) );
               SCIP_CALL( SCIPchgVarUbDive(scip, vars[i], SCIPround(scip, solval)) );

               SCIP_CALL( SCIPchgVarLbExactDive(scip, vars[i], newbound) );
               SCIP_CALL( SCIPchgVarUbExactDive(scip, vars[i], newbound) );

               SCIPrationalFreeBuffer(SCIPbuffer(scip), &newbound);
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
         SCIPwarningMessage(scip, "Error while solving LP in Exactsol Constraint Handler; exact LP solve terminated with code <%d>\n",retstat);
      }
#else
      SCIP_CALL( SCIPsolveExactDiveLP(scip, -1, &lperror, NULL) );
#endif

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
      /**@todo handle the unbounded case */
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &worksol) );
         conshdlrdata->nbufferedsols--;
      }
   }

   /* terminate exact diving */
   SCIP_CALL( SCIPendExactDive(scip) );

   if( foundsol )
   {
      SCIPdebugMessage("successfully found feasible improving solution, objval %g, upperbound %g\n",
         SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip));
      conshdlrdata->ncurrentstalls = 0;
   }
   else
   {
      SCIPdebugMessage("repaired solution not feasible or not improving, objval %g, upperbound %g \n",
         SCIPgetSolTransObj(scip, sol), SCIPgetUpperbound(scip));
      conshdlrdata->ncurrentstalls++;
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
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->hashedassignments, DEFAULT_SOLBUFSIZE) );
   SCIP_CALL( SCIPhashtableCreate(&(conshdlrdata->solhash), SCIPblkmem(scip), DEFAULT_SOLBUFSIZE, hashGetKeyAssignment, hashKeyEqAssignment, hashKeyValAssignment, NULL) );

   conshdlrdata->nhashedassignments = 0;
   conshdlrdata->lenhash = DEFAULT_SOLBUFSIZE;

   /* allocate data for solution buffer */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->solubuffer, DEFAULT_SOLBUFSIZE) );
   conshdlrdata->lensolubuffer = DEFAULT_SOLBUFSIZE;
   conshdlrdata->nbufferedsols = 0;

   conshdlrdata->ncurrentstalls = 0;
   conshdlrdata->probhasconteqs = -1;

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

   /* create exactsol constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );
   conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpExactSol, consEnfopsExactSol, consCheckExactSol, consLockExactSol,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* mark constraint handler as exact */
   SCIPconshdlrMarkExact(conshdlr);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyExactSol, NULL) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxExactSol) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitExactSol) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitExactSol) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeExactSol) );

   /* add exactsol constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/" CONSHDLR_NAME "/checkfpfeasibility",
         "should a solution be checked in floating-point arithmetic prior to being processed?",
         &conshdlrdata->checkfpfeasibility, TRUE, DEFAULT_CHECKFPFEASIBILITY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/maxstalls",
         "maximal number of consecutive repair calls without success",
         &conshdlrdata->maxstalls, TRUE, DEFAULT_MAXSTALLS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "constraints/" CONSHDLR_NAME "/solbufsize",
         "size of solution buffer",
         &conshdlrdata->solbufsize, TRUE, DEFAULT_SOLBUFSIZE, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "constraints/" CONSHDLR_NAME "/minimprove",
         "minimal percentage of primal improvement to trigger solution processing",
         &conshdlrdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
