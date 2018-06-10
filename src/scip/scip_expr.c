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

/**@file   scip_expr.c
 * @brief  public methods for expression handlers
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
#include "scip/visual.h"
#include "xml/xml.h"

#include "scip/scip_expr.h"
#include "scip/scip_mem.h"
#include "scip/scip_numerics.h"
#include "scip/scip_sol.h"
#include "scip/scip_var.h"

#include "scip/pub_message.h"
#include "scip/pub_nlp.h"
#include "scip/pub_var.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)

/** replaces array of variables in expression tree by corresponding transformed variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPgetExprtreeTransformedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetExprtreeTransformedVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPexprtreeGetNVars(tree) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetTransformedVars(scip, SCIPexprtreeGetNVars(tree), SCIPexprtreeGetVars(tree), SCIPexprtreeGetVars(tree)) );

   return SCIP_OKAY;
}

/** evaluates an expression tree for a primal solution or LP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPevalExprtreeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SOL*             sol,                /**< a solution, or NULL for current LP solution */
   SCIP_Real*            val                 /**< buffer to store value */
   )
{
   SCIP_Real* varvals;
   int nvars;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPevalExprtreeSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   nvars = SCIPexprtreeGetNVars(tree);

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPexprtreeEval(tree, NULL, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, SCIPexprtreeGetVars(tree), varvals) );

   SCIP_CALL( SCIPexprtreeEval(tree, varvals, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current global bounds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPevalExprtreeGlobalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
   )
{
   SCIP_INTERVAL* varvals;
   SCIP_VAR**     vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPevalExprtreeGlobalBounds", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   nvars = SCIPexprtreeGetNVars(tree);

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, NULL, val) );
      return SCIP_OKAY;
   }

   vars = SCIPexprtreeGetVars(tree);
   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -SCIPvarGetLbGlobal(vars[i])),  /*lint !e666*/
          infty2infty(SCIPinfinity(scip), infinity,  SCIPvarGetUbGlobal(vars[i]))); /*lint !e666*/
   }

   SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, varvals, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current local bounds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
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
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPevalExprtreeLocalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
   )
{
   SCIP_INTERVAL* varvals;
   SCIP_VAR**     vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPevalExprtreeLocalBounds", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   nvars = SCIPexprtreeGetNVars(tree);

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, NULL, val) );
      return SCIP_OKAY;
   }

   vars = SCIPexprtreeGetVars(tree);
   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      /* due to numerics, the lower bound on a variable in SCIP can be slightly higher than the upper bound
       * in this case, we take the most conservative way and switch the bounds
       * further, we translate SCIP's value for infinity to the users value for infinity
       */
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -MIN(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]))),  /*lint !e666*/
          infty2infty(SCIPinfinity(scip), infinity,  MAX(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])))); /*lint !e666*/
   }

   SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, varvals, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

#undef infty2infty
