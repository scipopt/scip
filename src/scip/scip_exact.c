/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   scip_exact.c
 * @brief  public methods for exact solving
 * @author Leon Eifler
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>
#ifndef _WIN32
#include <strings.h> /*lint --e{766}*/
#endif


#include "lpi/lpi.h"
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
#include "scip/lpexact_bounding.h"
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
#include "scip/sepastoreexact.h"
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

#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_exact.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"

#include "scip/pub_cons.h"
#include "scip/pub_fileio.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/pub_lpexact.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** enables or disables exact solving mode
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 */
SCIP_RETCODE SCIPenableExactSolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             enable              /**< enable exact solving (TRUE) or disable it (FALSE) */
   )
{
   assert(scip != NULL);

#ifndef SCIP_WITH_EXACTSOLVE
   if( enable )
   {
      SCIPerrorMessage("SCIP was compiled without exact solve support: cannot enable exact solving mode.\n");
      return SCIP_ERROR;
   }
#else
   /* skip if nothing has changed */
   if( enable == scip->set->exact_enable )
      return SCIP_OKAY;

   /* check stage and throw an error */
   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM )
   {
      SCIPerrorMessage("Exact solving mode can only be enabled/disabled before reading/creating a problem.\n");
      return SCIP_INVALIDCALL;
   }

   /* reoptimization in combination with exact solving has not been implemented */
   if( scip->set->reopt_enable )
   {
      SCIPerrorMessage("Exact solving mode not (yet) compatible with reoptimization.\n");
      return SCIP_INVALIDCALL;
   }

   scip->set->exact_enable = enable;
#endif

   return SCIP_OKAY;
}

/** returns whether the solution process should be probably correct
 *
 *  @return Returns TRUE if \SCIP is in exact solving mode, otherwise FALSE
 */
SCIP_Bool SCIPisExact(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (scip->set->exact_enable);
}

/** returns whether aggregation is allowed to use negative slack in exact solving mode
 *
 *  @return Returns TRUE if \SCIP is not in exact solving mode or aggregation is allowed to use negative slack
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_Bool SCIPallowNegSlack(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPallowNegSlack", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return (!SCIPisExact(scip)) || (scip->set->exact_allownegslack);
}

/** branches on an LP solution exactly; does not call branching rules, since fractionalities are assumed to small;
 *  if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPbranchLPExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPbranchLPExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecLPExact(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->origprob,
         scip->tree, scip->reopt, scip->lp, scip->branchcand, scip->eventqueue, scip->eventfilter,
         scip->primal->cutoffbound, TRUE, result) );

   return SCIP_OKAY;
}

/** adds row to exact separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        rowexact            /**< exact row to add */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddRowExact", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsepastoreExactAddCut(scip->sepastoreexact, scip->set, scip->eventqueue, rowexact) );

   return SCIP_OKAY;
}
