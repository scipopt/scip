/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   cons_integral.c
 * @ingroup DEFPLUGINS_CONS
 * @brief  constraint handler for the integrality constraint
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_integral.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/pub_sol.h"
#include "scip/scip_branch.h"
#include "scip/scip_cons.h"
#include "scip/scip_exact.h"
#include "scip/scip_lp.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_mem.h"
#include "scip/rational.h"
#include <string.h>


#define CONSHDLR_NAME          "integral"
#define CONSHDLR_DESC          "integrality constraint"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */


/** checks whether primal solution satisfies all integrality restrictions without tolerances */
static
SCIP_RETCODE checkIntegralityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to be checked */
   SCIP_Bool             printreason,        /**< should infeasisibility reason be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the lp enforcement call */
   )
{
   SCIP_RATIONAL* solval;
   SCIP_Bool integral;
   SCIP_VAR** vars;
   int nintegers;
   int ncontimplvars;
   int ncontvars;
   int v;

   assert(result != NULL);

   SCIPdebugMessage("checking exact integrality of primal solution:\n");

   integral = TRUE;

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &solval) );

   /* get all problem variables and integer region in vars array */
   SCIP_CALL( SCIPgetSolVarsData(scip, sol, &vars, &nintegers, NULL, NULL, NULL, NULL, &ncontimplvars, &ncontvars) );
   nintegers -= ncontimplvars + ncontvars;
   assert(nintegers >= 0);

   /* check whether primal solution satisfies all integrality restrictions */
   for( v = 0; v < nintegers && integral; ++v )
   {
      /* if the solution is exact we check the exact data, otherwise we check the fp data */
      if( SCIPsolIsExact(sol) )
         SCIPgetSolValExact(scip, sol, vars[v], solval);
      else
         SCIPrationalSetReal(solval, SCIPgetSolVal(scip, sol, vars[v]));

      assert(SCIPvarGetProbindex(vars[v]) == v);
      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY || SCIPvarGetType(vars[v]) == SCIP_VARTYPE_INTEGER );

      if( !SCIPrationalIsIntegral(solval) )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIPinfoMessage(scip, NULL, "violation: integrality condition of variable <%s> =",
               SCIPvarGetName(vars[v]));
            SCIPrationalMessage(SCIPgetMessagehdlr(scip), NULL, solval);
            SCIPinfoMessage(scip, NULL, "\n");
         }
         integral = FALSE;
      }
   }

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &solval);

   return SCIP_OKAY;
}

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

   SCIPdebugMsg(scip, "Enfolp method of integrality constraint: %d fractional variables\n", SCIPgetNLPBranchCands(scip));

    *result = SCIP_DIDNOTRUN;

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

   if( SCIPisExact(scip) && *result == SCIP_DIDNOTRUN )
   {
      SCIP_CALL( SCIPbranchLPExact(scip, result) );
   }

   /* if no branching was done, the LP solution was not fractional */
   if( *result == SCIP_DIDNOTRUN )
      *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxIntegral)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   int nbinvars;
   int nintvars;
   int i;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(conss == NULL);
   assert(nconss == 0);
   assert(result != NULL);

   SCIPdebugMsg(scip, "Enforelax method of integrality constraint\n");

   *result = SCIP_FEASIBLE;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, NULL, &nbinvars, &nintvars, NULL, NULL) );

   int nintegers = nbinvars + nintvars + SCIPgetNBinImplVars(scip) + SCIPgetNIntImplVars(scip);

   for( i = 0; i < nintegers; ++i )
   {
      assert(vars[i] != NULL);
      assert(SCIPvarIsIntegral(vars[i]));

      if( !SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[i])) )
      {
         if( SCIPisFeasEQ(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
         {
            SCIPdebugMsg(scip, "Cutoff for integral variable %s with bounds [%f, %f] and value %f\n", SCIPvarGetName(vars[i]),
                  SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]), SCIPgetSolVal(scip, sol, vars[i]));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else
         {
            /* @todo better way to handle this would be a BRANCHEXECRELAX callback that could also implement pseudo costs for
             * relaxation solutions instead of using the enforelaxcallback which is mainly intended for spatial branching
             */
            SCIP_CALL( SCIPaddExternBranchCand(scip, vars[i], 0.2, SCIPgetSolVal(scip, sol, vars[i])) );
            *result = SCIP_INFEASIBLE;
         }
      }
   }

   /* if we have found a branching candidate, immediately branch to be able to return SCIP_BRANCHED and stop the
    * enforcement loop
    */
   if( *result == SCIP_INFEASIBLE )
   {
      /* call branching methods for external candidates */
      SCIP_CALL( SCIPbranchExtern(scip, result) );

      /* since we only call it if we added external candidates, the branching rule should always be able to branch */
      assert(*result != SCIP_DIDNOTRUN);
   }

   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckIntegral)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real solval;
   int nintegers;
   int ncontimplvars;
   int ncontvars;
   int v;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIPdebugMsg(scip, "Check method of integrality constraint (checkintegrality=%u)\n", checkintegrality);

   *result = SCIP_FEASIBLE;

   if( !checkintegrality )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetSolVarsData(scip, sol, &vars, &nintegers, NULL, NULL, NULL, NULL, &ncontimplvars, &ncontvars) );
   nintegers -= ncontimplvars + ncontvars;
   assert(nintegers >= 0);

   if( !SCIPisExact(scip) )
   {
      for( v = 0; v < nintegers; ++v )
      {
         solval = SCIPgetSolVal(scip, sol, vars[v]);

         if( sol != NULL )
            SCIPupdateSolIntegralityViolation(scip, sol, EPSFRAC(solval, SCIPfeastol(scip)));

         if( !SCIPisFeasIntegral(scip, solval) )
         {
            *result = SCIP_INFEASIBLE;

            if( printreason )
            {
               SCIPinfoMessage(scip, NULL, "violation: integrality condition of variable <%s> = %.15g\n",
                  SCIPvarGetName(vars[v]), solval);
            }
            if( !completely )
               break;
         }
      }
   }
   /* in exact solving mode, we have to check integrality without tolerances */
   else
   {
      SCIP_CALL( checkIntegralityExact(scip, sol, printreason, result) );  
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockIntegral)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsIntegral)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real solval;
   SCIP_Real score;
   SCIP_Real bestscore;
   SCIP_Bool bestroundup;
   int nintegers;
   int ncontvars;
   int nbinimplvars;
   int nintimplvars;
   int ncontimplvars;
   int bestcandidx;
   int v;

   assert(scip != NULL);
   assert(diveset != NULL);
   assert(sol != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIPdebugMsg(scip, "integral constraint handler: determine diving bound changes\n");

   SCIP_CALL( SCIPgetSolVarsData(scip, sol, &vars, &nintegers, NULL, NULL,
         &nbinimplvars, &nintimplvars, &ncontimplvars, &ncontvars) );
   nintegers = nintegers - nbinimplvars - nintimplvars - ncontimplvars - ncontvars;
   assert(nintegers >= 0);

   bestscore = SCIP_REAL_MIN;
   bestcandidx = -1;
   *success = FALSE;
   bestroundup = FALSE; /* only for lint */

   /* loop over solution values and get score of fractional variables */
   for( v = 0; v < nintegers; ++v )
   {
      solval = SCIPgetSolVal(scip, sol, vars[v]);

      /* skip variable if solution value disagrees with the local bounds */
      if( ! SCIPisFeasIntegral(scip, solval) && SCIPisGE(scip, solval, SCIPvarGetLbLocal(vars[v])) && SCIPisLE(scip, solval, SCIPvarGetUbLocal(vars[v])) )
      {
         SCIP_Bool roundup;

         SCIP_CALL( SCIPgetDivesetScore(scip, diveset, SCIP_DIVETYPE_INTEGRALITY, vars[v], solval,
               solval - SCIPfloor(scip, solval), &score, &roundup) );

         /* we search for candidates with maximum score */
         if( score > bestscore )
         {
            bestcandidx = v;
            bestscore = score;
            bestroundup = roundup;
            *success = TRUE;
         }
      }
   }

   assert(!(*success) || bestcandidx >= 0);

   if( *success )
   {
      solval = SCIPgetSolVal(scip, sol, vars[bestcandidx]);

      /* if we want to round up the best candidate, it is added as the preferred bound change */
      SCIP_CALL( SCIPaddDiveBoundChange(scip, vars[bestcandidx], SCIP_BRANCHDIR_UPWARDS,
            SCIPceil(scip, solval), bestroundup) );
      SCIP_CALL( SCIPaddDiveBoundChange(scip, vars[bestcandidx], SCIP_BRANCHDIR_DOWNWARDS,
            SCIPfloor(scip, solval), ! bestroundup) );
   }

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
   SCIP_CONSHDLR* conshdlr;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpIntegral, consEnfopsIntegral, consCheckIntegral, consLockIntegral, NULL) );

   assert(conshdlr != NULL);

   /* mark constraint handler as exact */
   SCIPconshdlrMarkExact(conshdlr);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyIntegral, consCopyIntegral) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsIntegral) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxIntegral) );

   return SCIP_OKAY;
}
