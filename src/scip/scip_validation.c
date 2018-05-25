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

/**@file   scip_validation.c
 * @brief  public methods for validation
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

#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_validation.h"

#include "scip/pub_message.h"
#include "scip/pub_misc.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** validate the result of the solve
 *
 *  the validation includes
 *
 *  - checking the feasibility of the incumbent solution in the original problem (using SCIPcheckSolOrig())
 *
 *  - checking if the objective bounds computed by SCIP agree with external primal and dual reference bounds.
 *
 *  All external reference bounds the original problem space and the original objective sense.
 *
 *  For infeasible problems, +/-SCIPinfinity() should be passed as reference bounds depending on the objective sense
 *  of the original problem.
 */
SCIP_RETCODE SCIPvalidateSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             primalreference,    /**< external primal reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             dualreference,      /**< external dual reference value for the problem, or SCIP_UNKNOWN */
   SCIP_Real             reftol,             /**< relative tolerance for acceptable violation of reference values */
   SCIP_Bool             quiet,              /**< TRUE if no status line should be printed */
   SCIP_Bool*            feasible,           /**< pointer to store if the best solution is feasible in the original problem,
                                               *  or NULL */
   SCIP_Bool*            primalboundcheck,   /**< pointer to store if the primal bound respects the given dual reference
                                               *  value, or NULL */
   SCIP_Bool*            dualboundcheck      /**< pointer to store if the dual bound respects the given primal reference
                                               *  value, or NULL */
   )
{
   SCIP_Bool localfeasible;
   SCIP_Bool localprimalboundcheck;
   SCIP_Bool localdualboundcheck;
   SCIP_Real primviol;
   SCIP_Real dualviol;
   assert(scip != NULL);

   /* if no problem exists, there is no need for validation */
   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      if( feasible != NULL )
         *feasible = TRUE;
      if( primalboundcheck != NULL )
         *primalboundcheck = TRUE;
      if( dualboundcheck != NULL )
         *dualboundcheck = TRUE;

      return SCIP_OKAY;
   }

   localfeasible = TRUE;
   localdualboundcheck = TRUE;

   /* check the best solution for feasibility in the original problem */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_SOL* bestsol = SCIPgetBestSol(scip);
      SCIP_Real checkfeastolfac;
      SCIP_Real oldfeastol;

      assert(bestsol != NULL);

      /* scale feasibility tolerance by set->num_checkfeastolfac */
      oldfeastol = SCIPfeastol(scip);
      SCIP_CALL( SCIPgetRealParam(scip, "numerics/checkfeastolfac", &checkfeastolfac) );
      if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
      {
         SCIP_CALL( SCIPchgFeastol(scip, oldfeastol * checkfeastolfac) );
      }

      SCIP_CALL( SCIPcheckSolOrig(scip, bestsol, &localfeasible, !quiet, TRUE) );

      /* restore old feasibilty tolerance */
      if( !SCIPisEQ(scip, checkfeastolfac, 1.0) )
      {
         SCIP_CALL( SCIPchgFeastol(scip, oldfeastol) );
      }
   }
   else
   {
      localfeasible = TRUE;
   }

   primviol = 0.0;
   dualviol = 0.0;
   /* check the primal and dual bounds computed by SCIP against the external reference values within reference tolerance */
   /* solution for an infeasible problem */
   if( SCIPgetNSols(scip) > 0 && ((SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE && SCIPisInfinity(scip, dualreference))
            || (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE && SCIPisInfinity(scip, -dualreference))) )
      localprimalboundcheck = FALSE;
   else
   {
      /* check if reference primal bound is not better than the proven dual bound and, if SCIP claims to be optimal,
       * if the
       */
      SCIP_Real pb = SCIPgetPrimalbound(scip);
      SCIP_Real db = SCIPgetDualbound(scip);

      /* compute the relative violation between the primal bound and dual reference value, and vice versa */
      if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
      {
         if( dualreference != SCIP_UNKNOWN ) /*lint !e777 */
            primviol = SCIPrelDiff(dualreference, pb);
         if( primalreference != SCIP_UNKNOWN ) /*lint !e777 */
            dualviol = SCIPrelDiff(db, primalreference);
      }
      else
      {
         if( dualreference != SCIP_UNKNOWN ) /*lint !e777 */
            primviol = SCIPrelDiff(pb, dualreference);

         if( primalreference != SCIP_UNKNOWN ) /*lint !e777 */
            dualviol = SCIPrelDiff(primalreference, db);
      }
      primviol = MAX(primviol, 0.0);
      dualviol = MAX(dualviol, 0.0);

      localprimalboundcheck = EPSP(reftol, primviol);
      localdualboundcheck = EPSP(reftol, dualviol);
   }

   if( !quiet )
   {
      SCIPinfoMessage(scip, NULL, "Validation         : ");
      if( ! localfeasible )
         SCIPinfoMessage(scip, NULL, "Fail (infeasible)");
      else if( ! localprimalboundcheck )
         SCIPinfoMessage(scip, NULL, "Fail (primal bound)");
      else if( ! localdualboundcheck )
         SCIPinfoMessage(scip, NULL, "Fail (dual bound)");
      else
         SCIPinfoMessage(scip, NULL, "Success");
      SCIPinfoMessage(scip, NULL, "\n");
      SCIPinfoMessage(scip, NULL, "  %-17s: %10u\n", "cons violation", !localfeasible);
      SCIPinfoMessage(scip, NULL, "  %-17s: %10.8g (reference: %16.9e)\n", "primal violation", primviol, dualreference);
      SCIPinfoMessage(scip, NULL, "  %-17s: %10.8g (reference: %16.9e)\n", "dual violation", dualviol, primalreference);
   }

   if( feasible != NULL )
      *feasible = localfeasible;
   if( primalboundcheck != NULL )
      *primalboundcheck = localprimalboundcheck;
   if( dualboundcheck != NULL )
      *dualboundcheck = localdualboundcheck;

   return SCIP_OKAY;
}
