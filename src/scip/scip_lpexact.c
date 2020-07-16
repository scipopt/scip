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

/**@file   scip_lpexact.c
 * @brief  public methods for the exact LP relaxation, rows and columns
 * @author Leon Eifler
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


#include "lpiexact/lpiexact.h"
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
#include "scip/lpexact.h"
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
#include "scip/rational.h"
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
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"

#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_tree.h"
#include "scip/scip_message.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif


/** decreases usage counter of LP row, and frees memory if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPreleaseRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT**       row                 /**< pointer to LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPreleaseRowExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIProwExactRelease(row, scip->mem->probmem, scip->set, scip->lpexact) );

   return SCIP_OKAY;
}

static
void SCIProwExactForceSort(
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SET*             set                 /**< SCIP settings */
   )
{
   assert(row != NULL);
   assert(row->nunlinked == row->len);
   assert(row->nlpcols == 0);

   SCIPsetDebugMsg(set, "merging row <%s>\n", row->fprow->name);

   /* do nothing on empty rows; if row is sorted, nothing has to be done */
   if( row->len > 0 && (!row->lpcolssorted || !row->nonlpcolssorted) )
   {
      SCIP_COLEXACT** cols;
      int* cols_index;
      SCIP_Rational** vals;
      int s;
      int t;

      /* make sure, the row is sorted */
      SCIProwExactSort(row);
      assert(row->lpcolssorted);
      assert(row->nonlpcolssorted);

      /* merge equal columns, thereby recalculating whether the row's activity is always integral */
      cols = row->cols;
      cols_index = row->cols_index;
      vals = row->vals;
      assert(cols != NULL);
      assert(cols_index != NULL);
      assert(vals != NULL);

      t = 0;
      row->integral = TRUE;
      assert(!RatIsZero(vals[0]));
      assert(row->linkpos[0] == -1);

      for( s = 1; s < row->len; ++s )
      {
         assert(!RatIsZero(vals[s]));
         assert(row->linkpos[s] == -1);

         if( cols[s] == cols[t] )
         {
            /* merge entries with equal column */
            RatAdd(vals[t], vals[t], vals[s]);
         }
         else
         {
            /* go to the next entry, overwriting current entry if coefficient is zero */
            if( !RatIsZero(vals[t]) )
            {
               row->integral = row->integral && SCIPcolIsIntegral(cols[t]->fpcol) && RatIsIntegral(vals[t]);
               t++;
            }
            cols[t] = cols[s];
            cols_index[t] = cols_index[s];
            RatSet(vals[t], vals[s]);
         }
      }
      if( !RatIsZero(vals[t]) )
      {
         row->integral = row->integral && SCIPcolIsIntegral(cols[t]->fpcol) && RatIsIntegral(vals[t]);
         t++;
      }
      assert(s == row->len);
      assert(t <= row->len);

      row->len = t;
      row->nunlinked = t;
   }

#ifndef NDEBUG
   /* check for double entries */
   {
      int i;
      int j;

      for( i = 0; i < row->len; ++i )
      {
         assert(row->cols[i] != NULL);
         assert(row->cols[i]->index == row->cols_index[i]);
         for( j = i+1; j < row->len; ++j )
            assert(row->cols[i] != row->cols[j]);
      }
   }
#endif
}

/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @attention If a coefficients absolute value is below the SCIP epsilon tolerance, the variable with its value is not added.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarsToRowEx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Rational**       vals                /**< values of coefficients */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarsToRowEx", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   SCIP_CALL( SCIProwExactEnsureSize(row, scip->mem->probmem, scip->set, SCIProwGetNNonz(row->fprow) + nvars) );

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarAddToRowExact(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lpexact,
            row, vals[v]) );
   }
   SCIPdebug(SCIProwPrint(row->fprow, SCIPgetMessagehdlr(scip), NULL));
   SCIPdebug(SCIProwExactPrint(row, SCIPgetMessagehdlr(scip), NULL));

   /* force the row sorting */
   SCIProwExactForceSort(row, scip->set);
   row->fprow->rowexact = row;

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients from a constraint handler
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateEmptyRowConsExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT**       rowexact,           /**< pointer to row */
   SCIP_ROW*             fprow,              /**< corresponding fp-row */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs                 /**< right hand side of row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowConsExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreateExact(rowexact, fprow, scip->mem->probmem, scip->set,
                                 scip->stat, scip->lpexact, 0, NULL, NULL, lhs, rhs,
                                 SCIP_ROWORIGINTYPE_CONS, SCIProwGetOriginCons(fprow)) );

   return SCIP_OKAY;
}

/** returns the feasibility of a row for the given primal solution
 *
 *  @return the feasibility of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void SCIPgetRowSolFeasibilityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowSolFeasibilityExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      SCIProwExactGetSolFeasibility(row, scip->set, scip->stat, sol, result);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      SCIProwExactGetLPFeasibility(row, scip->set, scip->stat, scip->lpexact, result);
   else
      SCIProwExactGetPseudoFeasibility(row, scip->set, scip->stat, result);
}

/** returns the activity of a row for the given primal solution
 *
 *  @return the activitiy of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
void SCIPgetRowSolActivityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             useexact,           /**< true if sol should be considered instead of sol */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowSolActivityExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      SCIProwExactGetSolActivity(row, scip->set, scip->stat, sol, useexact, result);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      RatSet(result, SCIProwExactGetLPActivity(row, scip->set, scip->stat, scip->lpexact));
   else
      RatSet(result, SCIProwExactGetPseudoActivity(row, scip->set, scip->stat));
}


/** output exact row to file stream via the message handler system
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPprintRowex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintRowex", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIProwExactPrint(row, scip->messagehdlr, file);

   return SCIP_OKAY;
}

/** gets objective value of current exact LP (which is the sum of column and loose objective value)
 *
 *  @return the objective value of current LP (which is the sum of column and loose objective value).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This method returns the objective value of the current LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the LP solution status returned by
 *        SCIPgetLPSolstat() is SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Rational* SCIPgetLPExactObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Rational* res;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPlpExactGetObjval(scip->lpexact, scip->set, scip->transprob, res);


   return res;
}


/** returns whether the exact lp was solved */
SCIP_Bool SCIPlpExactIsSolved(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->lpexact != NULL);

   return scip->lpexact->solved && scip->lpexact->flushed;
}


/** initiates LP diving, making methods SCIPchgVarObjDive(), SCIPchgVarLbDive(), and SCIPchgVarUbDive() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note diving is allowed even if the current LP is not flushed, not solved, or not solved to optimality; be aware
 *  that solving the (first) diving LP may take longer than expect and that the latter two cases could stem from
 *  numerical troubles during the last LP solve; because of this, most users will want to call this method only if
 *  SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL
 */
SCIP_RETCODE SCIPstartExactDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstartDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   assert(SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_FOCUSNODE);

   if( SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("already in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot start exact diving while being in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   // TODO: maybe replace with an exact version?
   if( !SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      SCIPerrorMessage("cannot start diving if LP has not been constructed\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeHasCurrentNodeLP(scip->tree));

   SCIP_CALL( SCIPlpExactStartDive(scip->lpexact, scip->mem->probmem, scip->set, scip->stat) );

   /* TODO: I don't think we need that in exact diving mode, but we have to check! */
   /* remember the relaxation solution to reset it later */
   //if( SCIPisRelaxSolValid(scip) )
   //{
   //   SCIP_CALL( SCIPtreeStoreRelaxSol(scip->tree, scip->set, scip->relaxation, scip->transprob) );
   //}

   return SCIP_OKAY;
}

/** quits exact LP diving and resets bounds and objective values of columns to the current node's values
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPendExactDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPendDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("not in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* unmark the diving flag in the LP and reset all variables' objective and bound values */
   SCIP_CALL( SCIPlpExactEndDive(scip->lpexact, scip->mem->probmem, scip->set, scip->messagehdlr, scip->stat, scip->eventqueue, scip->eventfilter,
         scip->transprob, scip->transprob->vars, scip->transprob->nvars) );

   /* the lower bound may have changed slightly due to LP resolve in SCIPlpEndDive() */
   if( !scip->lpexact->resolvelperror && scip->tree->focusnode != NULL && SCIPlpIsRelax(scip->lp) && SCIPlpExactIsSolved(scip) )
   {
      assert(SCIPtreeIsFocusNodeLPConstructed(scip->tree));
      SCIP_CALL( SCIPnodeUpdateLowerboundLP(scip->tree->focusnode, scip->set, scip->stat, scip->tree, scip->transprob,
            scip->origprob, scip->lp) );
   }
   /* reset the probably changed LP's cutoff bound */
   SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->transprob, scip->primal->cutoffbound) );
   assert(scip->lp->cutoffbound == scip->primal->cutoffbound); /*lint !e777*/

   /* if a new best solution was created, the cutoff of the tree was delayed due to diving;
    * the cutoff has to be done now.
    */
   if( scip->tree->cutoffdelayed )
   {
      SCIP_CALL( SCIPtreeCutoff(scip->tree, scip->reopt, scip->mem->probmem, scip->set, scip->stat, scip->eventfilter,
            scip->eventqueue, scip->lp, scip->primal->cutoffbound) );
   }

   /* TODO: see comment in SCIPstartExactDive */
   /* if a relaxation was stored before diving, restore it now */
   //if( scip->tree->probdiverelaxstored )
   //{
   //   SCIP_CALL( SCIPtreeRestoreRelaxSol(scip->tree, scip->set, scip->relaxation, scip->transprob) );
   //}

   return SCIP_OKAY;
}

/** solves the exact LP of the current dive; no separation or pricing is applied
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note be aware that the LP solve may take longer than expected if SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL,
 *  compare the explanation of SCIPstartDive()
 */
SCIP_RETCODE SCIPsolveExactDiveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the diving LP was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveDiveLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lpexact) )
   {
      SCIPerrorMessage("not in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( cutoff != NULL )
      *cutoff = FALSE;

   /* solve diving LP */
   SCIP_CALL( SCIPlpExactSolveAndEval(scip->lpexact, scip->lp, scip->set, scip->messagehdlr, scip->mem->probmem, scip->stat,
         scip->eventqueue, scip->eventfilter, scip->transprob, (SCIP_Longint)itlim, lperror, FALSE) );

   /* the LP is infeasible or the objective limit was reached */
   if( SCIPlpExactGetSolstat(scip->lpexact) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpExactGetSolstat(scip->lpexact) == SCIP_LPSOLSTAT_OBJLIMIT
      || (SCIPlpExactGetSolstat(scip->lpexact) == SCIP_LPSOLSTAT_OPTIMAL &&
         RatIsGE(SCIPgetLPExactObjval(scip), SCIPgetCutoffboundExact(scip))) )
   {
      /* analyze the infeasible LP (only if the objective was not changed, all columns are in the LP, and no external
       * pricers exist) */
      //if( !(SCIPlpDivingObjChanged(scip->lp) || SCIPlpDivingRowsChanged(scip->lp))
      //   && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) )
      //{
      //   SCIP_CALL( SCIPconflictAnalyzeLP(scip->conflict, scip->conflictstore, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
      //         scip->origprob, scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->cliquetable, NULL) );
      //}

      if( cutoff != NULL )
         *cutoff = TRUE;
   }

   return SCIP_OKAY;
}

/** changes variable's lower bound in current exact dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgVarLbExactDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Rational*        newbound            /**< new value for bound */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("not in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPvarChgLbExactDive(var, scip->set, scip->lpexact, newbound) );

   return SCIP_OKAY;
}

/** changes variable's upper bound in current exact dive
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPchgVarUbExactDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Rational*        newbound            /**< new value for bound */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("not in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPvarChgUbExactDive(var, scip->set, scip->lpexact, newbound) );

   return SCIP_OKAY;
}
