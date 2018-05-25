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

/**@file   scip_nlp.c
 * @brief  public methods for nonlinear relaxations
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

#include "scip/scip_mem.h"
#include "scip/scip_nlp.h"
#include "scip/scip_sol.h"

#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_nlp.h"
#include "scip/pub_paramset.h"


/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif

/** method to call, when the priority of an NLPI was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNlpiPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSetPriorityNlpi() to mark the nlpis unsorted */
   SCIP_CALL( SCIPsetNlpiPriority(scip, (SCIP_NLPI*)paramdata, SCIPparamGetInt(param)) );

   return SCIP_OKAY;
}
/** includes an NLPI in SCIP */
SCIP_RETCODE SCIPincludeNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi                /**< NLPI data structure */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(nlpi != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPincludeNlpi", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether NLPI is already present */
   if( SCIPfindNlpi(scip, SCIPnlpiGetName(nlpi)) != NULL )
   {
      SCIPerrorMessage("NLPI <%s> already included.\n", SCIPnlpiGetName(nlpi));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPsetIncludeNlpi(scip->set, nlpi) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "nlpi/%s/priority", SCIPnlpiGetName(nlpi));
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of NLPI <%s>", SCIPnlpiGetName(nlpi));
   SCIP_CALL( SCIPaddIntParam(scip, paramname, paramdesc,
         NULL, FALSE, SCIPnlpiGetPriority(nlpi), INT_MIN/4, INT_MAX/4,
         paramChgdNlpiPriority, (SCIP_PARAMDATA*)nlpi) ); /*lint !e740*/

   /* pass message handler (may be NULL) */
   SCIP_CALL( SCIPnlpiSetMessageHdlr(nlpi, scip->messagehdlr) );

   return SCIP_OKAY;
}

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_NLPI* SCIPfindNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of NLPI */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(name != NULL);

   return SCIPsetFindNlpi(scip->set, name);
}

/** returns the array of currently available NLPIs (sorted by priority) */
SCIP_NLPI** SCIPgetNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSortNlpis(scip->set);

   return scip->set->nlpis;
}

/** returns the number of currently available NLPIs */
int SCIPgetNNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nnlpis;
}

/** sets the priority of an NLPI */
SCIP_RETCODE SCIPsetNlpiPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of the NLPI */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPsetSetPriorityNlpi(scip->set, nlpi, priority);

   return SCIP_OKAY;
}

/** returns whether the NLP relaxation has been enabled
 *
 *  If the NLP relaxation is enabled, then SCIP will construct the NLP relaxation when the solving process is about to begin.
 *  To check whether an NLP is existing, use SCIPisNLPConstructed().
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @see SCIPenableNLP
 */
SCIP_Bool SCIPisNLPEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisNLPEnabled", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->transprob->nlpenabled;
}

/** marks that there are constraints that are representable by nonlinear rows
 *
 *  This method should be called by a constraint handler if it has constraints that have a representation as nonlinear rows.
 *
 *  The function should be called before the branch-and-bound process is initialized, e.g., when presolve is exiting.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
void SCIPenableNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPenableNLP", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   scip->transprob->nlpenabled = TRUE;
}

/** returns, whether an NLP has been constructed
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Bool SCIPisNLPConstructed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisNLPConstructed", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return (scip->nlp != NULL);
}

/** returns whether the NLP has a continuous variable in a nonlinear term
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Bool SCIPhasNLPContinuousNonlinearity(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPhasNLPContinuousNonlinearity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }

   return SCIPnlpHasContinuousNonlinearity(scip->nlp);
}

/** gets current NLP variables along with the current number of NLP variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store the array of NLP variables, or NULL */
   int*                  nvars               /**< pointer to store the number of NLP variables, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPVarsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      if( vars != NULL )
         *vars = SCIPnlpGetVars(scip->nlp);
      if( nvars != NULL )
         *nvars = SCIPnlpGetNVars(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   return SCIP_OKAY;
}

/** gets array with variables of the NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_VAR** SCIPgetNLPVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPVars", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   return SCIPnlpGetVars(scip->nlp);
}

/** gets current number of variables in NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNNLPVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNLPVars", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }

   return SCIPnlpGetNVars(scip->nlp);
}

/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPVarsNonlinearity(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcount             /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPVarsNonlinearity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpGetVarsNonlinearity(scip->nlp, nlcount) );

   return SCIP_OKAY;
}

/** returns dual solution values associated with lower bounds of NLP variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real* SCIPgetNLPVarsLbDualsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPVarsLbDualsol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   return SCIPnlpGetVarsLbDualsol(scip->nlp);
}

/** returns dual solution values associated with upper bounds of NLP variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real* SCIPgetNLPVarsUbDualsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPVarsUbDualsol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   return SCIPnlpGetVarsUbDualsol(scip->nlp);
}

/** gets current NLP nonlinear rows along with the current number of NLP nonlinear rows
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPNlRowsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW***         nlrows,             /**< pointer to store the array of NLP nonlinear rows, or NULL */
   int*                  nnlrows             /**< pointer to store the number of NLP nonlinear rows, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPNlRowsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   if( nlrows != NULL )
      *nlrows = SCIPnlpGetNlRows(scip->nlp);
   if( nnlrows != NULL )
      *nnlrows = SCIPnlpGetNNlRows(scip->nlp);

   return SCIP_OKAY;
}

/** gets array with nonlinear rows of the NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NLROW** SCIPgetNLPNlRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPNlRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return NULL; /*lint !e527*/
   }

   return SCIPnlpGetNlRows(scip->nlp);
}

/** gets current number of nonlinear rows in NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNNLPNlRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNLPNlRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return 0; /*lint !e527*/
   }

   return SCIPnlpGetNNlRows(scip->nlp);
}

/** adds a nonlinear row to the NLP. This row is captured by the NLP.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< nonlinear row to add to NLP */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpAddNlRow(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, nlrow) );

   return SCIP_OKAY;
}

/** makes sure that the NLP of the current node is flushed
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPflushNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPflushNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpFlush(scip->nlp, scip->mem->probmem, scip->set) );

   return SCIP_OKAY;
}

/** sets or clears initial primal guess for NLP solution (start point for NLP solver)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNLPInitialGuess(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            initialguess        /**< values of initial guess (corresponding to variables from SCIPgetNLPVarsData), or NULL to use no start point */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNLPInitialGuess", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpSetInitialGuess(scip->nlp, SCIPblkmem(scip), initialguess) );

   return SCIP_OKAY;
}

/** sets initial primal guess for NLP solution (start point for NLP solver)
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNLPInitialGuessSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution which values should be taken as initial guess, or NULL for LP solution */
   )
{
   SCIP_Real* vals;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNLPInitialGuessSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, SCIPnlpGetNVars(scip->nlp)) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPnlpGetNVars(scip->nlp), SCIPnlpGetVars(scip->nlp), vals) );
   SCIP_CALL( SCIPnlpSetInitialGuess(scip->nlp, SCIPblkmem(scip), vals) );
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

/** solves the current NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsolveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpSolve(scip->nlp, SCIPblkmem(scip), scip->set, scip->messagehdlr, scip->stat) );

   return SCIP_OKAY;
}

/** gets solution status of current NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NLPSOLSTAT SCIPgetNLPSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPSolstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return SCIP_NLPSOLSTAT_UNKNOWN; /*lint !e527*/
   }

   return SCIPnlpGetSolstat(scip->nlp);
}

/** gets termination status of last NLP solve
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_NLPTERMSTAT SCIPgetNLPTermstat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPTermstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return SCIP_NLPTERMSTAT_OTHER; /*lint !e527*/
   }

   return SCIPnlpGetTermstat(scip->nlp);
}

/** gives statistics (number of iterations, solving time, ...) of last NLP solve
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpGetStatistics(scip->nlp, statistics) );

   return SCIP_OKAY;
}

/** gets objective value of current NLP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Real SCIPgetNLPObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetObjval(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALID;
   }
}

/** indicates whether a feasible solution for the current NLP is available
 * thus, returns whether the solution status <= feasible
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_Bool SCIPhasNLPSolution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPhasNLPSolution", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }

   return SCIPnlpHasSolution(scip->nlp);
}

/** gets fractional variables of last NLP solution along with solution values and fractionalities
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPFracVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           fracvars,           /**< pointer to store the array of NLP fractional variables, or NULL */
   SCIP_Real**           fracvarssol,        /**< pointer to store the array of NLP fractional variables solution values, or NULL */
   SCIP_Real**           fracvarsfrac,       /**< pointer to store the array of NLP fractional variables fractionalities, or NULL */
   int*                  nfracvars,          /**< pointer to store the number of NLP fractional variables , or NULL */
   int*                  npriofracvars       /**< pointer to store the number of NLP fractional variables with maximal branching priority, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPFracVars", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpGetFracVars(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, fracvars, fracvarssol, fracvarsfrac, nfracvars, npriofracvars) );

   return SCIP_OKAY;
}

/** gets integer parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPIntPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPIntPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpGetIntPar(scip->nlp, type, ival) );

   return SCIP_OKAY;
}

/** sets integer parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNLPIntPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNLPIntPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpSetIntPar(scip->nlp, type, ival) );

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPRealPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPRealPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpGetRealPar(scip->nlp, type, dval) );

   return SCIP_OKAY;
}

/** sets floating point parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNLPRealPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNLPRealPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpSetRealPar(scip->nlp, type, dval) );

   return SCIP_OKAY;
}

/** gets string parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPStringPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPStringPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpGetStringPar(scip->nlp, type, sval) );

   return SCIP_OKAY;
}

/** sets string parameter of NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNLPStringPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNLPStringPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpSetStringPar(scip->nlp, type, sval) );

   return SCIP_OKAY;
}

/** writes current NLP to a file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPwriteNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpWrite(scip->nlp, scip->set, scip->messagehdlr, filename) );

   return SCIP_OKAY;
}

/** gets the NLP interface and problem used by the SCIP NLP;
 *  with the NLPI and its problem you can use all of the methods defined in nlpi/nlpi.h;
 *
 *  @warning You have to make sure, that the full internal state of the NLPI does not change or is recovered completely
 *           after the end of the method that uses the NLPI. In particular, if you manipulate the NLP or its solution
 *           (e.g. by calling one of the SCIPnlpiAdd...() or the SCIPnlpiSolve() method), you have to check in advance
 *           whether the NLP is currently solved.  If this is the case, you have to make sure, the internal solution
 *           status is recovered completely at the end of your method. Additionally you have to resolve the NLP with
 *           SCIPnlpiSolve() in order to reinstall the internal solution status.
 *
 *  @warning Make also sure, that all parameter values that you have changed are set back to their original values.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNLPI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI**           nlpi,               /**< pointer to store the NLP solver interface */
   SCIP_NLPIPROBLEM**    nlpiproblem         /**< pointer to store the NLP solver interface problem */
   )
{
   assert(nlpi != NULL);
   assert(nlpiproblem != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNLPI", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   *nlpi = SCIPnlpGetNLPI(scip->nlp);
   *nlpiproblem = SCIPnlpGetNLPIProblem(scip->nlp);

   return SCIP_OKAY;
}


/*
 * NLP diving methods
 */

/**@name NLP Diving Methods */
/**@{ */

/** initiates NLP diving making methods SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), SCIPchgVarsBoundsDiveNLP(), and SCIPsolveDiveNLP() available
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPstartDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPstartDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpStartDive(scip->nlp, SCIPblkmem(scip), scip->set) );

   return SCIP_OKAY;
}

/** ends NLP diving
 *
 *  Resets changes made by SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), and SCIPchgVarsBoundsDiveNLP().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPendDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPendDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpEndDive(scip->nlp, SCIPblkmem(scip), scip->set) );

   return SCIP_OKAY;
}

/** changes linear objective coefficient of a variable in diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarObjDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new value for coefficient */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarObjDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpChgVarObjDive(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, var, coef) );

   return SCIP_OKAY;
}

/** changes bounds of a variable in diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarBoundsDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which bounds to change */
   SCIP_Real             lb,                 /**< new lower bound */
   SCIP_Real             ub                  /**< new upper bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarBoundsDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpChgVarBoundsDive(scip->nlp, var, lb, ub) );

   return SCIP_OKAY;
}

/** changes bounds of a set of variables in diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgVarsBoundsDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables which bounds to changes */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds */
   SCIP_Real*            ubs                 /**< new upper bounds */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarsBoundsDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpChgVarsBoundsDive(scip->nlp, scip->set, nvars, vars, lbs, ubs) );

   return SCIP_OKAY;
}

/** solves diving NLP
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsolveDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsolveDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlpSolveDive(scip->nlp, SCIPblkmem(scip), scip->set, scip->messagehdlr, scip->stat) );

   return SCIP_OKAY;
}

/**@} */


/*
 * NLP nonlinear row methods
 */

/** creates and captures an NLP row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   int                   nquadvars,          /**< number variables in quadratic terms */
   SCIP_VAR**            quadvars,           /**< variables in quadratic terms, or NULL if nquadvars == 0 */
   int                   nquadelems,         /**< number of elements in quadratic term */
   SCIP_QUADELEM*        quadelems,          /**< elements (i.e., monomials) in quadratic term, or NULL if nquadelems == 0 */
   SCIP_EXPRTREE*        expression,         /**< nonlinear expression, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_EXPRCURV         curvature           /**< curvature of the nonlinear row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowCreate(nlrow, scip->mem->probmem, scip->set,
         name, constant, nlinvars, linvars, lincoefs, nquadvars, quadvars, nquadelems, quadelems, expression, lhs, rhs, curvature) );

   return SCIP_OKAY;
}

/** creates and captures an NLP nonlinear row without any coefficients
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateEmptyNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< pointer to nonlinear row */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateEmptyNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowCreate(nlrow, scip->mem->probmem, scip->set,
         name, 0.0, 0, NULL, NULL, 0, NULL, 0, NULL, NULL, lhs, rhs, SCIP_EXPRCURV_UNKNOWN) );

   return SCIP_OKAY;
}

/** creates and captures an NLP row from a linear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcreateNlRowFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< pointer to nonlinear row */
   SCIP_ROW*             row                 /**< the linear row to copy */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateNlRowFromRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowCreateFromRow(nlrow, scip->mem->probmem, scip->set, row) );

   return SCIP_OKAY;
}

/** increases usage counter of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPcaptureNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcaptureNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIPnlrowCapture(nlrow);

   return SCIP_OKAY;
}

/** decreases usage counter of NLP nonlinear row, and frees memory if necessary
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_RETCODE SCIPreleaseNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow               /**< pointer to nonlinear row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPreleaseNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowRelease(nlrow, scip->mem->probmem, scip->set) );

   return SCIP_OKAY;
}

/** changes left hand side of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgNlRowLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgNlRowLhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgLhs(nlrow, scip->set, scip->stat, scip->nlp, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgNlRowRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgNlRowRhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgRhs(nlrow, scip->set, scip->stat, scip->nlp, rhs) );

   return SCIP_OKAY;
}

/** changes constant of NLP nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgNlRowConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real             constant            /**< new value for constant */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgNlRowConstant", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgConstant(nlrow, scip->set, scip->stat, scip->nlp, constant) );

   return SCIP_OKAY;
}

/** adds variable with a linear coefficient to the nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddLinearCoefToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value of coefficient in linear part of row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddLinearCoefToNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowAddLinearCoef(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, var, val) );

   return SCIP_OKAY;
}

/** adds variables with linear coefficients to the row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddLinearCoefsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real*            vals                /**< values of coefficients in linear part of row */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddLinearCoefsToNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPnlrowAddLinearCoef(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, vars[v], vals[v]) );
   }

   return SCIP_OKAY;
}

/** changes linear coefficient of a variables in a row
 *
 *  Setting the coefficient to 0.0 means that it is removed from the row
 *  the variable does not need to exists before.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgNlRowLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   assert(var != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgNlRowLinearCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgLinearCoef(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, var, coef) );

   return SCIP_OKAY;
}

/** adds quadratic variable to the nonlinear row
 *
 *  After adding a quadratic variable, it can be used to add quadratic elements.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddQuadVarToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddQuadVarToNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, scip->mem->probmem, scip->set, var) );

   return SCIP_OKAY;
}

/** adds quadratic variables to the nonlinear row
 *
 *  After adding quadratic variables, they can be used to add quadratic elements.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddQuadVarsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars                /**< problem variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddQuadVarsToNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowEnsureQuadVarsSize(nlrow, scip->mem->probmem, scip->set, SCIPnlrowGetNQuadVars(nlrow) + nvars) );
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, scip->mem->probmem, scip->set, vars[v]) );
   }

   return SCIP_OKAY;
}

/** add a quadratic element to the nonlinear row
 *
 *  Variable indices of the quadratic element need to be relative to quadratic variables array of row.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddQuadElementToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_QUADELEM         quadelem            /**< quadratic element */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddQuadElementToNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, quadelem) );

   /* invalidate curvature */
   if( quadelem.coef != 0.0 )
      SCIPnlrowSetCurvature(nlrow, SCIP_EXPRCURV_UNKNOWN);

   return SCIP_OKAY;
}

/** adds quadratic elements to the nonlinear row
 *
 *  Variable indices of the quadratic elements need to be relative to quadratic variables array of row.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddQuadElementsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements */
   )
{
   int v;

   assert(nquadelems == 0 || quadelems != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddQuadElementsToNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowEnsureQuadElementsSize(nlrow, scip->mem->probmem, scip->set, SCIPnlrowGetNQuadElems(nlrow) + nquadelems) );
   for( v = 0; v < nquadelems; ++v )
   {
      SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, quadelems[v]) );
   }

   /* invalidate curvature */
   SCIPnlrowSetCurvature(nlrow, SCIP_EXPRCURV_UNKNOWN);

   return SCIP_OKAY;
}

/** changes coefficient in quadratic part of a row
 *
 *  Setting the coefficient in the quadelement to 0.0 means that it is removed from the row
 *  the element does not need to exists before.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPchgNlRowQuadElement(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_QUADELEM         quadelement         /**< new quadratic element, or update for existing one */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgNlRowQuadElement", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgQuadElem(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, quadelement) );

   return SCIP_OKAY;
}

/** sets or deletes expression tree in the nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNlRowExprtree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_EXPRTREE*        exprtree            /**< expression tree, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNlRowExprtree", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgExprtree(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, exprtree) );

   /* invalidate curvature */
   SCIPnlrowSetCurvature(nlrow, SCIP_EXPRCURV_UNKNOWN);

   return SCIP_OKAY;
}

/** sets a parameter of expression tree in the nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNlRowExprtreeParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   paramidx,           /**< index of parameter in expression tree */
   SCIP_Real             paramval            /**< new value of parameter in expression tree */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNlRowExprtreeParam", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgExprtreeParam(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, paramidx, paramval) );

   return SCIP_OKAY;
}

/** sets parameters of expression tree in the nonlinear row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPsetNlRowExprtreeParams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real*            paramvals           /**< new values of parameter in expression tree */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPsetNlRowExprtreeParams", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgExprtreeParams(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, paramvals) );

   return SCIP_OKAY;
}

/** recalculates the activity of a nonlinear row in the last NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPrecalcNlRowNLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrecalcNlRowNLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("do not have NLP for computing NLP activity\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlrowRecalcNLPActivity(nlrow, scip->set, scip->stat, scip->nlp) );

   return SCIP_OKAY;
}

/** returns the activity of a nonlinear row in the last NLP solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowNLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            activity            /**< pointer to store activity value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowNLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("do not have NLP for computing NLP activity\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, scip->set, scip->stat, scip->nlp, activity) );

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row in the last NLP solution: negative value means infeasibility
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowNLPFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            feasibility         /**< pointer to store feasibility value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowNLPFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("do not have NLP for computing NLP feasibility\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, scip->set, scip->stat, scip->nlp, feasibility) );

   return SCIP_OKAY;
}

/** recalculates the activity of a nonlinear row for the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPrecalcNlRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrecalcNlRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowRecalcPseudoActivity(nlrow, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** gives the activity of a nonlinear row for the current pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            pseudoactivity      /**< pointer to store pseudo activity value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, scip->set, scip->stat, pseudoactivity) );

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row for the current pseudo solution: negative value means infeasibility
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowPseudoFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            pseudofeasibility   /**< pointer to store pseudo feasibility value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowPseudoFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowGetPseudoFeasibility(nlrow, scip->set, scip->stat, pseudofeasibility) );

   return SCIP_OKAY;
}

/** recalculates the activity of a nonlinear row in the last NLP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPrecalcNlRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPrecalcNlRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowRecalcNLPActivity(nlrow, scip->set, scip->stat, scip->nlp) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowRecalcPseudoActivity(nlrow, scip->set, scip->stat) );
   }

   return SCIP_OKAY;
}

/** gives the activity of a nonlinear row in the last NLP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            activity            /**< pointer to store activity value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, scip->set, scip->stat, scip->nlp, activity) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, scip->set, scip->stat, activity) );
   }

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row in the last NLP or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            feasibility         /**< pointer to store feasibility value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, scip->set, scip->stat, scip->nlp, feasibility) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoFeasibility(nlrow, scip->set, scip->stat, feasibility) );
   }

   return SCIP_OKAY;
}

/** gives the activity of a nonlinear row for the given primal solution or NLP solution or pseudo solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowSolActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for NLP solution of pseudo solution */
   SCIP_Real*            activity            /**< pointer to store activity value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowSolActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
   {
      SCIP_CALL( SCIPnlrowGetSolActivity(nlrow, scip->set, scip->stat, sol, activity) );
   }
   else if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, scip->set, scip->stat, scip->nlp, activity) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, scip->set, scip->stat, activity) );
   }

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row for the given primal solution
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowSolFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< pointer to store feasibility value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
   {
      SCIP_CALL( SCIPnlrowGetSolFeasibility(nlrow, scip->set, scip->stat, sol, feasibility) );
   }
   else if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, scip->set, scip->stat, scip->nlp, feasibility) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoFeasibility(nlrow, scip->set, scip->stat, feasibility) );
   }

   return SCIP_OKAY;
}

/** gives the minimal and maximal activity of a nonlinear row w.r.t. the variable's bounds
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPgetNlRowActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPgetNlRowActivityBounds", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowGetActivityBounds(nlrow, scip->set, scip->stat, minactivity, maxactivity) );

   return SCIP_OKAY;
}

/** output nonlinear row to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPprintNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(nlrow != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowPrint(nlrow, scip->messagehdlr, file) );

   return SCIP_OKAY;
}
