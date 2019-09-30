
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_lpex.c
 * @brief  public methods for the LP relaxation, rows and columns
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


#include "lpi/lpiex.h"
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
#include "scip/varex.h"
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
   SCIP_ROWEX**          row                 /**< pointer to LP row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPreleaseRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIProwexRelease(row, scip->mem->probmem, scip->set, scip->lpex) );

   return SCIP_OKAY;
}

static
void SCIProwexForceSort(
   SCIP_ROWEX*           row,
   SCIP_SET*             set
   )
{
   assert(row != NULL);
   assert(row->nunlinked == row->len);
   assert(row->nlpcols == 0);

   SCIPsetDebugMsg(set, "merging row <%s>\n", row->name);

   /* do nothing on empty rows; if row is sorted, nothing has to be done */
   if( row->len > 0 && (!row->lpcolssorted || !row->nonlpcolssorted) )
   {
      SCIP_COLEX** cols;
      int* cols_index;
      SCIP_Rational** vals;
      int s;
      int t;

      /* make sure, the row is sorted */
      SCIProwexSort(row);
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
      assert(!RisZero(vals[0]));
      assert(row->linkpos[0] == -1);

      for( s = 1; s < row->len; ++s )
      {
         assert(!RisZero(vals[s]));
         assert(row->linkpos[s] == -1);

         if( cols[s] == cols[t] )
         {
            /* merge entries with equal column */
            Radd(vals[t], vals[t], vals[s]);
         }
         else
         {
            /* go to the next entry, overwriting current entry if coefficient is zero */
            if( !RisZero(vals[t]) )
            {
               row->integral = row->integral && SCIPcolIsIntegral(cols[t]->fpcol) && RisIntegral(vals[t]);
               t++;
            }
            cols[t] = cols[s];
            cols_index[t] = cols_index[s];
            Rset(vals[t], vals[s]);
         }
      }
      if( !RisZero(vals[t]) )
      {
         row->integral = row->integral && SCIPcolIsIntegral(cols[t]->fpcol) && RisIntegral(vals[t]);
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
   SCIP_ROWEX*           row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Rational**       vals                /**< values of coefficients */
   )
{
   SCIP_Real* realvals;
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddVarsToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   //SCIP_CALL( SCIPallocBufferArray(scip, &realvals, nvars) );
   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   SCIP_CALL( SCIProwexEnsureSize(row, scip->mem->probmem, scip->set, SCIProwGetNNonz(row->fprow) + nvars) );

   /* todo exip: do we need the sorting? */
   /* delay the row sorting */
   //SCIProwDelaySort(row);

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarAddToRowExact(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lpex,
            row, vals[v]) );
      //realvals[v] = RgetRealApprox(vals[v]);
   }
   SCIPdebug(SCIProwPrint(row->fprow, SCIPgetMessagehdlr(scip), NULL));
   SCIPdebug(SCIProwexPrint(row, SCIPgetMessagehdlr(scip), NULL));

   /* force the row sorting */
   SCIProwexForceSort(row, scip->set);
   SCIPhashtableInsert(scip->lpex->exrowhash, (void*) row);

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
   SCIP_ROWEX**          rowex,              /**< pointer to row */
   SCIP_ROW*             fprow,              /**< corresponding fp-row */
   SCIP_Rational*        lhs,                /**< left hand side of row */
   SCIP_Rational*        rhs                 /**< right hand side of row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyRowCons", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreateExact(rowex, fprow, scip->mem->probmem, scip->set,
                                 scip->stat, scip->lpex, 0, NULL, NULL, lhs, rhs,
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
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      SCIProwexGetSolFeasibility(row, scip->set, scip->stat, sol, result);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      SCIProwexGetLPFeasibility(row, scip->set, scip->stat, scip->lpex, result);
   else
      SCIProwexGetPseudoFeasibility(row, scip->set, scip->stat, result);
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
   SCIP_ROWEX*           row,                /**< LP row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_SOLEX*           solex,              /**< primal CIP solution */
   SCIP_Bool             useexact,           /**< true if solex should be considered instead of sol */
   SCIP_Rational*        result              /**< result pointer */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetRowSolActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      SCIProwexGetSolActivity(row, scip->set, scip->stat, sol, solex, useexact, result);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      Rset(result, SCIProwexGetLPActivity(row, scip->set, scip->stat, scip->lpex));
   else
      Rset(result, SCIProwexGetPseudoActivity(row, scip->set, scip->stat));
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
   SCIP_ROWEX*           row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIProwexPrint(row, scip->messagehdlr, file);

   return SCIP_OKAY;
}


SCIP_Bool SCIPlpexIsSolved(
   SCIP* scip
   )
{
   assert(scip != NULL);
   assert(scip->lpex != NULL);

   return scip->lpex->solved && scip->lpex->flushed;
}