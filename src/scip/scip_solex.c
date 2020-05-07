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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_sol.c
 * @brief  public methods for solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Robert Lion Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif


#include "lpi/lpi.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"
#include "scip/benders.h"
#include "scip/benderscut.h"
#include "scip/branch.h"
#include "scip/branch_nodereopt.h"
#include "scip/clock.h"
#include "scip/compr.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/cons_linear.h"
#include "scip/cutpool.h"
#include "scip/cuts.h"
#include "scip/debug.h"
#include "scip/def.h"
#include "scip/dialog.h"
#include "scip/dialog_default.h"
#include "scip/disp.h"
#include "scip/event.h"
#include "scip/heur.h"
#include "scip/heur_ofins.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heuristics.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/interrupt.h"
#include "scip/lp.h"
#include "scip/lpex.h"
#include "scip/mem.h"
#include "scip/message_default.h"
#include "scip/misc.h"
#include "scip/nlp.h"
#include "scip/nodesel.h"
#include "scip/paramset.h"
#include "scip/presol.h"
#include "scip/presolve.h"
#include "scip/pricer.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/prop.h"
#include "scip/reader.h"
#include "scip/relax.h"
#include "scip/reopt.h"
#include "scip/retcode.h"
#include "scip/scipbuildflags.h"
#include "scip/scipcoreplugins.h"
#include "scip/scipgithash.h"
#include "scip/sepa.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/syncstore.h"
#include "scip/table.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/varex.h"
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solex.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"

#include "scip/pub_cons.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/pub_varex.h"

#include "scip/struct_sol.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** checks solution for exact feasibility in original problem without adding it to the solution store; to improve the
 *  performance we use the following order when checking for violations:
 *
 *  1. variable bounds
 *  2. constraint handlers with positive or zero priority that don't need constraints (e.g. integral constraint handler)
 *  3. original constraints
 *  4. constraint handlers with negative priority that don't need constraints (e.g. Benders' decomposition constraint handler)
 */
static
SCIP_RETCODE checkSolexOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< Should the reason for the violation be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool             checkmodifiable     /**< have modifiable constraint to be checked? */
   )
{
   SCIP_Rational* solval;
   SCIP_Rational* lb;
   SCIP_Rational* ub;
   SCIP_RESULT result;
   int v;
   int c;
   int h;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(feasible != NULL);
   assert(SCIPisExactSol(scip, sol));

   SCIP_CALL( SCIPcheckStage(scip, "checkSolexOrig", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   *feasible = TRUE;

   SCIPsolResetViolations(sol);

   if( !printreason )
      completely = FALSE;

   RatCreateBuffer(SCIPbuffer(scip), &solval);

   /* check bounds */
   if( checkbounds )
   {
      for( v = 0; v < scip->origprob->nvars; ++v )
      {
         SCIP_VAR* var;

         var = scip->origprob->vars[v];
         SCIPsolGetValExact(solval, sol, scip->set, scip->stat, var);

         lb = SCIPvarGetLbOriginalExact(var);
         ub = SCIPvarGetUbOriginalExact(var);

         if( RatIsLT(solval, lb) || RatIsGT(solval, ub) )
         {
            *feasible = FALSE;

            if( printreason )
            {
               SCIPmessagePrintInfo(scip->messagehdlr, "solution violates original bounds of variable <%s> [%g,%g] solution value <%g>\n",
                  SCIPvarGetName(var), RatApproxReal(lb), RatApproxReal(ub), RatApproxReal(solval));
            }

            if( !completely )
            {
               RatFreeBuffer(SCIPbuffer(scip), &solval);
               return SCIP_OKAY;
            }
         }
      }
   }

   RatFreeBuffer(SCIPbuffer(scip), &solval);

   /* call constraint handlers with positive or zero check priority that don't need constraints */
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      if( SCIPconshdlrGetCheckPriority(scip->set->conshdlrs[h]) >= 0 )
      {
         if( !SCIPconshdlrNeedsCons(scip->set->conshdlrs[h]) )
         {
            SCIP_CALL( SCIPconshdlrCheck(scip->set->conshdlrs[h], scip->mem->probmem, scip->set, scip->stat, sol,
                  checkintegrality, checklprows, printreason, completely, &result) );

            if( result != SCIP_FEASIBLE )
            {
               *feasible = FALSE;

               if( !completely )
                  return SCIP_OKAY;
            }
         }
      }
      /* constraint handlers are sorted by priority, so we can break when reaching the first one with negative priority */
      else
         break;
   }

   /* check original constraints
    *
    * in general modifiable constraints can not be checked, because the variables to fulfill them might be missing in
    * the original problem; however, if the solution comes from a heuristic during presolving modifiable constraints
    * have to be checked;
    */
   for( c = 0; c < scip->origprob->nconss; ++c )
   {
      if( SCIPconsIsChecked(scip->origprob->conss[c]) && (checkmodifiable || !SCIPconsIsModifiable(scip->origprob->conss[c])) )
      {
         /* check solution */
         SCIP_CALL( SCIPconsCheck(scip->origprob->conss[c], scip->set, sol,
               checkintegrality, checklprows, printreason, &result) );

         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;

            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   /* call constraint handlers with negative check priority that don't need constraints;
    * continue with the first constraint handler with negative priority which caused us to break in the above loop */
   for( ; h < scip->set->nconshdlrs; ++h )
   {
      assert(SCIPconshdlrGetCheckPriority(scip->set->conshdlrs[h]) < 0);
      if( !SCIPconshdlrNeedsCons(scip->set->conshdlrs[h]) )
      {
         SCIP_CALL( SCIPconshdlrCheck(scip->set->conshdlrs[h], scip->mem->probmem, scip->set, scip->stat, sol,
               checkintegrality, checklprows, printreason, completely, &result) );

         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;

            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution
 *
 *  @return value of variable in primal CIP solution, or in current LP/pseudo solution
 *
 *  @pre In case the solution pointer @p sol is @b NULL, that means it is asked for the LP or pseudo solution, this method
 *       can only be called if @p scip is in the solving stage \ref SCIP_STAGE_SOLVING. In any other case, this method
 *       can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
void SCIPgetSolexVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_VAR*             var,                /**< variable to get value for */
   SCIP_Rational*        res                 /**< resulting rational */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolexVal", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert( var->scip == scip );

   if( sol != NULL )
   {
      SCIPsolGetValExact(res, sol, scip->set, scip->stat, var);
   }
   else
   {
      SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolexVal(sol==NULL)", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

      RatSet(res, SCIPvarGetSolex(var, SCIPtreeHasCurrentNodeLP(scip->tree)));
   }
}

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @return transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 */
void SCIPgetSolexTransObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo objective value */
   SCIP_Rational*        res                 /**< result pointer to store rational */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolexTransObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( sol != NULL )
      RatSet(res, SCIPsolGetObjExact(sol, scip->set, scip->transprob, scip->origprob));
   else
   {
      SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetSolexTransObj(sol==NULL)", \
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
      if( SCIPtreeHasCurrentNodeLP(scip->tree) )
         SCIPlpexGetObjval(scip->lpex, scip->set, scip->transprob, res);
      else
         SCIPlpexGetPseudoObjval(scip->lpex, scip->set, scip->transprob, res);
   }
}

/** overwrite the fp-values in a solution with the rounded exact ones */
SCIP_RETCODE SCIPoverwriteFPsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   int v;
   int nvars;
   SCIP_Rational* res;
   SCIP_VAR** vars;

   assert(scip != NULL);
   assert(sol != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPoverwriteFPsol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIP_CALL( SCIPsolOverwriteFPSolExact(sol, scip->set, scip->stat, scip->origprob, scip->transprob, scip->tree) );

   return SCIP_OKAY;
}

/** print a rational solution */
SCIP_RETCODE SCIPprintSolex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Rational* objvalue;
   SCIP_Bool currentsol;
   SCIP_Bool oldquiet = FALSE;

   assert(SCIPisTransformed(scip) || sol != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintSolex", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   currentsol = (sol == NULL);
   if( currentsol )
   {
      SCIP_CALL( SCIPcheckStage(scip, "SCIPprintSolex(sol==NULL)", \
            FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

      /* create a temporary solution that is linked to the current solution */
      SCIP_CALL( SCIPsolCreateCurrentSolExact(&sol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
            scip->tree, scip->lpex, NULL) );
   }

   if( file != NULL && scip->messagehdlr != NULL )
   {
      oldquiet = SCIPmessagehdlrIsQuiet(scip->messagehdlr);
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, FALSE);
   }

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "objective value:                 ");

   /** @todo exip: convert to origobj, or extern objval respectively */
   if( SCIPsolIsOriginal(sol) )
      objvalue = SCIPsolGetOrigObjExact(sol);
   else
      objvalue = SCIPsolGetObjExact(sol, scip->set, scip->transprob, scip->origprob);

   RatMessage(scip->messagehdlr, file, objvalue);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");

   SCIP_CALL( SCIPsolPrintExact(sol, scip->set, scip->messagehdlr, scip->stat, scip->origprob, scip->transprob, file, FALSE,
         printzeros) );

   if( file != NULL && scip->messagehdlr != NULL )
   {
      SCIPmessagehdlrSetQuiet(scip->messagehdlr, oldquiet);
   }

   if( currentsol )
   {
      /* free temporary solution */
      SCIP_CALL( SCIPsolFree(&sol, scip->mem->probmem, scip->primal) );
   }

   return SCIP_OKAY;
}

/** returns TRUE if the solution is an exact rational solution */
SCIP_Bool SCIPisExactSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(SCIPisTransformed(scip) || sol != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPisExactSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsolIsExact(sol);
}

/** checks solution for exact feasibility without adding it to the solution store
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcheckSolex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< Should all reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked if printreason is true? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            feasible            /**< stores whether given solution is feasible */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcheckSolex", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* return immediately if the solution is of type partial */
   if( SCIPsolIsPartial(sol) )
   {
      SCIPerrorMessage("Cannot check feasibility of partial solutions.");
      return SCIP_INVALIDDATA;
   }

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || scip->set->misc_exactsolve;

   if( !printreason )
      completely = FALSE;

   if( SCIPsolIsOriginal(sol) )
   {
      /* SCIPsolCheck() can only be called on transformed solutions */
      SCIP_CALL( checkSolexOrig(scip, sol, feasible, printreason, completely, checkbounds, checkintegrality, checklprows, FALSE) );
   }
   else
   {
      SCIP_CALL( SCIPsolCheck(sol, scip->set, scip->messagehdlr, scip->mem->probmem, scip->stat, scip->transprob,
            printreason, completely, checkbounds, checkintegrality, checklprows, feasible) );
   }

   return SCIP_OKAY;
}

/** checks exact solution for feasibility in original problem without adding it to the solution store;
 *  this method is used to double check a solution in order to validate the presolving process
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcheckSolexOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool             completely          /**< Should all violations be checked if printreason is true? */
   )
{
   assert(scip != NULL);
   assert(sol != NULL);
   assert(feasible != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcheckSolexOrig", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* return immediately if the solution is of type partial */
   if( SCIPsolIsPartial(sol) )
   {
      SCIPerrorMessage("Cannot check feasibility of partial solutions.");
      return SCIP_INVALIDDATA;
   }

   if( !printreason )
      completely = FALSE;

   /* check solution in original problem; that includes bounds, integrality, and non modifiable constraints */
   SCIP_CALL( checkSolexOrig(scip, sol, feasible, printreason, completely, TRUE, TRUE, TRUE, FALSE) );

   return SCIP_OKAY;
}