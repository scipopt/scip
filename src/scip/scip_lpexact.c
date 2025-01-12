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
#include "scip/exprinterpret.h"
#include "scip/nlpi.h"
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
#include "scip/sepastoreexact.h"
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
#include "scip/scip_lpexact.h"
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

/** increases usage counter of exact LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcaptureRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row                 /**< row to capture */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcaptureRowExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIProwExactCapture(row);

   return SCIP_OKAY;
}

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

/** changes left hand side of exact LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgRowExactLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        lhs                 /**< new left hand side */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgRowExactLhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(!SCIPlpExactDiving(scip->lpexact) || (row->lppos == -1));

   SCIP_CALL( SCIProwExactChgLhs(row, scip->set, scip->lpexact, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of exact LP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgRowExactRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_Rational*        rhs                 /**< new right hand side */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgRowExactRhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(!SCIPlpExactDiving(scip->lpexact) || (row->lppos == -1));

   SCIP_CALL( SCIProwExactChgRhs(row, scip->set, scip->lpexact, rhs) );

   return SCIP_OKAY;
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
SCIP_RETCODE SCIPaddVarsToRowExact(
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

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarsToRowExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   SCIP_CALL( SCIProwExactEnsureSize(row, scip->mem->probmem, scip->set, SCIProwGetNNonz(row->fprow) + nvars) );

   /* delay the row sorting */
   SCIProwExactDelaySort(row);

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
   SCIP_ROW*             fprowrhs,           /**< rhs-part of fp-relaxation of this row if necessary, NULL otherwise */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs,                /**< right hand side of row */
   SCIP_Bool             isfprelaxable       /**< is it possible to make fp-relaxation of this row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowConsExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwExactCreate(rowexact, fprow, fprowrhs, scip->mem->probmem, scip->set,
                                 scip->stat, scip->lpexact, 0, NULL, NULL, lhs, rhs,
                                 SCIP_ROWORIGINTYPE_CONS, isfprelaxable, SCIProwGetOriginCons(fprow)) );

   return SCIP_OKAY;
}

/** creates and captures an exact LP row without any coefficients from a separator
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateEmptyRowExactSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT**       rowexact,           /**< pointer to exact row */
   SCIP_ROW*             fprow,              /**< corresponding fp approximation/relaxation */
   SCIP_SEPA*            sepa,               /**< separator that creates the row */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs,                /**< right hand side of row */
   SCIP_Bool             hasfprelaxation     /**< the the fprow a relaxation or only an approximation of the exact row? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowSepa", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwExactCreate(rowexact, fprow, NULL, scip->mem->probmem, scip->set, scip->stat,
         scip->lpexact, 0, NULL, NULL, lhs, rhs, SCIP_ROWORIGINTYPE_SEPA, hasfprelaxation, (void*) sepa) );

   return SCIP_OKAY;
}

/** creates and captures an exact LP row from an existing fp row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateRowExactFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             fprow               /**< corresponding fp approximation/relaxation */
   )
{
   assert(fprow != NULL);
   assert(fprow->rowexact == NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateRowExactFromRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIProwExactCreateFromRow(fprow, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lpexact) );

   return SCIP_OKAY;
}

/** generates two fprows that are a relaxation of the exact row wrt the lhs/rhs, respectively
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgenerateFpRowsFromRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< SCIP exact row */
   SCIP_ROW*             rowlhs,             /**< fp row-relaxation wrt lhs */
   SCIP_ROW*             rowrhs,             /**< fp row-relaxation wrt rhs */
   SCIP_Bool*            onerowrelax,        /**< is one row enough to represent the exact row */
   SCIP_Bool*            hasfprelax          /**< is it possible to generate relaxations at all for this row? */
   )
{
   assert(row != NULL);
   assert(rowlhs != NULL);
   assert(rowrhs != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgenerateFpRowsFromRowExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIProwExactGenerateFpRows(scip->mem->probmem, scip->set, scip->stat,  scip->eventqueue, scip->lpexact, scip->transprob, row, rowlhs, rowrhs, onerowrelax, hasfprelax) );

   return SCIP_OKAY;
}

/** returns the feasibility of a row for the given primal solution
 *
 *  @return the feasibility of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetRowSolFeasibilityExact(
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
      SCIP_CALL( SCIProwExactGetLPFeasibility(row, scip->set, scip->stat, scip->lpexact, result) );
   else
      SCIP_CALL( SCIProwExactGetPseudoFeasibility(row, scip->set, scip->stat, result) );

   return SCIP_OKAY;
}

/** returns the activity of a row for the given primal solution
 *
 *  @return the activitiy of a row for the given primal solution
 *
 *  @pre this method can be called in one of the following stages of the SCIP solving process:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetRowSolActivityExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             useexact,           /**< true if sol should be considered instead of sol */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowSolActivityExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      SCIP_CALL( SCIProwExactGetSolActivity(row, scip->set, scip->stat, sol, useexact, result) );
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      RatSet(result, SCIProwExactGetLPActivity(row, scip->stat, scip->lpexact));
   else
      RatSet(result, SCIProwExactGetPseudoActivity(row, scip->stat));

   return SCIP_OKAY;
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
SCIP_RETCODE SCIPprintRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintRowExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIProwExactPrint(row, scip->messagehdlr, file);

   return SCIP_OKAY;
}

/** gets objective value of current exact LP (which is the sum of column and loose objective value)
 *
 *  @return the objective value of current exact LP (which is the sum of column and loose objective value).
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note This method returns the objective value of the current exact LP solution, which might be primal or dual infeasible
 *        if a limit was hit during solving. It must not be used as a dual bound if the exact LP solution status returned by
 *        SCIPgetLPExactSolstat() is SCIP_LPSOLSTAT_ITERLIMIT or SCIP_LPSOLSTAT_TIMELIMIT.
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
void SCIPgetLPExactObjval(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPExactObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPlpExactGetObjval(scip->lpexact, scip->set, result);
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

/** gets solution status of current exact LP
  *
  *  @return the solution status of current exact LP.
  *
  *  @pre This method can be called if @p scip is in one of the following stages:
  *       - \ref SCIP_STAGE_SOLVING
  *
  *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
  */
SCIP_LPSOLSTAT SCIPgetLPExactSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLPExactSolstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* We can check the floating point flag here since the exact and floating point LP is constructed at the same
    * time.
    */
   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
       return SCIPlpExactGetSolstat(scip->lpexact);
   else
       return SCIP_LPSOLSTAT_NOTSOLVED;
}

/** initiates exact LP diving, making methods SCIPchgVarObjExactDive(), SCIPchgVarLbExactDive(), and SCIPchgVarUbExactDive() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note In parallel to exact LP diving, this method also starts the regular LP diving mode by calling SCIPstartDive().
 */
SCIP_RETCODE SCIPstartExactDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPstartExactDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   assert(SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_FOCUSNODE);

   if( SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("already in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("cannot start exact diving while being in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot start exact diving while being in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* We start the exact LP dive parallel to the floating-point LP dive. This is necessary because we need to work with
    * several flags and counters of the floating-point LP.
    */
   SCIP_CALL( SCIPstartDive(scip) );

   SCIP_CALL( SCIPlpExactStartDive(scip->lpexact, scip->mem->probmem, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** checks if exact diving mode is possible at this point in time
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note In parallel to exact LP diving, this method also starts the regular LP diving mode by calling SCIPstartDive().
 */
SCIP_Bool SCIPisExactDivePossible(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisExactDivePossible", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   assert(SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_FOCUSNODE);

   if( SCIPlpExactDiving(scip->lpexact) )
      return FALSE;

   if( SCIPlpDiving(scip->lp) )
      return FALSE;

   if( SCIPtreeProbing(scip->tree) )
      return FALSE;

   if( !SCIPlpIsSolved(scip->lp) )
      return FALSE;

   return TRUE;
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

   SCIP_CALL( SCIPcheckStage(scip, "SCIPendExactDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("not in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /** @todo exip: adress problem when user calls `SCIPendDive` in between */
   /* end floating-point LP dive, see comment in SCIPstartExactDive() */
   SCIP_CALL( SCIPendDive(scip) );

   /* unmark the diving flag in the exact LP and reset all variables' objective and bound values */
   SCIP_CALL( SCIPlpExactEndDive(scip->lpexact, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
         scip->transprob->vars, scip->transprob->nvars) );

   /* reset the probably changed LP's cutoff bound */
   SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->transprob, scip->primal->cutoffbound) );
   assert(scip->lp->cutoffbound == scip->primal->cutoffbound); /*lint !e777*/

   /* we have to set the exact diving flag temporarilly to TRUE since SCIPendDive() needs to know that this happend
    * in exact diving mode */
   scip->lpexact->diving = TRUE;

   scip->lpexact->diving = FALSE;

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
 */
SCIP_RETCODE SCIPsolveExactDiveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror,            /**< pointer to store whether an unresolved LP error occurred */
   SCIP_Bool*            cutoff              /**< pointer to store whether the diving LP was infeasible or the objective
                                              *   limit was reached (or NULL, if not needed) */
   )
{
   SCIP_Rational* objval;

   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveExactDiveLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("not in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( cutoff != NULL )
      *cutoff = FALSE;

   /* solve diving LP */
   SCIP_CALL( SCIPlpExactSolveAndEval(scip->lpexact, scip->lp, scip->set, scip->messagehdlr, scip->mem->probmem, scip->stat,
         scip->eventqueue, scip->transprob, (SCIP_Longint)itlim, lperror, FALSE) );

   if( !(*lperror) )
   {
      SCIP_CALL( RatCreateBuffer(scip->set->buffer, &objval) );
      SCIPgetLPExactObjval(scip, objval);

      /* the LP is infeasible or the objective limit was reached */
      if( SCIPlpExactGetSolstat(scip->lpexact) == SCIP_LPSOLSTAT_INFEASIBLE || SCIPlpExactGetSolstat(scip->lpexact) == SCIP_LPSOLSTAT_OBJLIMIT
         || (SCIPlpExactGetSolstat(scip->lpexact) == SCIP_LPSOLSTAT_OPTIMAL &&
            RatIsGE(objval, SCIPgetCutoffboundExact(scip))) )
      {
         if( cutoff != NULL )
            *cutoff = TRUE;
      }

      RatFreeBuffer(scip->set->buffer, &objval);
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

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarLbExactDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

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

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarUbExactDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPlpExactDiving(scip->lpexact) )
   {
      SCIPerrorMessage("not in exact diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPvarChgUbExactDive(var, scip->set, scip->lpexact, newbound) );

   return SCIP_OKAY;
}

/** writes current exact LP to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPwriteLPexact(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   SCIP_Bool cutoff;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteLPexact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      SCIP_CALL( SCIPconstructCurrentLP(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
            scip->tree, scip->reopt, scip->lp, scip->pricestore, scip->sepastore, scip->cutpool, scip->branchcand,
            scip->eventqueue, scip->eventfilter, scip->cliquetable, FALSE, &cutoff) );
   }

   /* we need a flushed lp to write the current lp */
   SCIP_CALL( SCIPsepastoreExactSyncLPs(scip->sepastoreexact, scip->mem->probmem, scip->set, scip->lpexact, scip->eventqueue) );
   SCIP_CALL( SCIPlpExactFlush(scip->lpexact, scip->mem->probmem, scip->set, scip->eventqueue) );

   SCIP_CALL( SCIPlpExactWrite(scip->lpexact, filename) );

   return SCIP_OKAY;
}
