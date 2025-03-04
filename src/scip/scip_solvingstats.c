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

/**@file   scip_solvingstats.c
 * @ingroup OTHER_CFILES
 * @brief  public methods for querying solving statistics
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Leona Gottwald
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Kati Wolter
 * @author Mohammed Ghannam
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/branch.h"
#include "scip/clock.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/conflict.h"
#include "scip/conflictstore.h"
#include "scip/debug.h"
#include "scip/disp.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/pricestore.h"
#include "scip/primal.h"
#include "scip/prob.h"
#include "scip/pub_benderscut.h"
#include "scip/pub_benders.h"
#include "scip/pub_branch.h"
#include "scip/pub_compr.h"
#include "scip/pub_cons.h"
#include "scip/pub_cutpool.h"
#include "scip/pub_cutsel.h"
#include "scip/pub_expr.h"
#include "scip/pub_heur.h"
#include "scip/pub_history.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_nlpi.h"
#include "scip/pub_presol.h"
#include "scip/pub_pricer.h"
#include "scip/pub_prop.h"
#include "scip/pub_reader.h"
#include "scip/pub_relax.h"
#include "scip/pub_reopt.h"
#include "scip/pub_sepa.h"
#include "scip/pub_sol.h"
#include "scip/pub_table.h"
#include "scip/pub_var.h"
#include "scip/reader.h"
#include "scip/reopt.h"
#include "scip/scip_benders.h"
#include "scip/scip_datatree.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_table.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"
#include "scip/sepastore.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/stat.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/syncstore.h"
#include "scip/table.h"
#include "scip/tree.h"
#include "scip/var.h"
#include <string.h>

/** gets number of branch and bound runs performed, including the current run
 *
 *  @return the number of branch and bound runs performed, including the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
int SCIPgetNRuns(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNRuns", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nruns;
}

/** gets number of reoptimization runs performed, including the current run
 *
 *  @return the number of reoptimization runs performed, including the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
int SCIPgetNReoptRuns(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNReoptRuns", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nreoptruns;
}

/** add given number to the number of processed nodes in current run and in all runs, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
void SCIPaddNNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Longint          nnodes              /**< number of processed nodes to add to the statistics */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPaddNNodes", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   scip->stat->nnodes += nnodes;
   scip->stat->ntotalnodes += nnodes;
}

/** gets number of processed nodes in current run, including the focus node
 *
 *  @return the number of processed nodes in current run, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_Longint SCIPgetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNodes", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nnodes;
}

/** gets total number of processed nodes in all runs, including the focus node
 *
 *  @return the total number of processed nodes in all runs, including the focus node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_Longint SCIPgetNTotalNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNTotalNodes", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->ntotalnodes;
}

/** gets number of leaf nodes processed with feasible relaxation solution
 *
 * @return number of leaf nodes processed with feasible relaxation solution
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_Longint SCIPgetNFeasibleLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNFeasibleLeaves", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nfeasleaves;
}

/** gets number of infeasible leaf nodes processed
 *
 * @return number of infeasible leaf nodes processed
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_Longint SCIPgetNInfeasibleLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNInfeasibleLeaves", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->ninfeasleaves;
}

/** gets number of processed leaf nodes that hit LP objective limit
 *
 * @return number of processed leaf nodes that hit LP objective limit
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_Longint SCIPgetNObjlimLeaves(
   SCIP*                 scip                /**< Scip data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNObjlimLeaves", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nobjleaves;
}

/** gets number of global bound changes
 *
 * @return number of global bound changes
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
int SCIPgetNRootboundChgs(
   SCIP*                 scip                /**< Scip data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNRootboundChgs", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nrootboundchgs;
}

/** gets number of global bound changes applied in the current run
 *
 * @return number of global bound changes
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
int SCIPgetNRootboundChgsRun(
   SCIP*                 scip                /**< Scip data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNRootboundChgsRun", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nrootboundchgsrun;
}

/** gets number of times a selected node was from a cut off subtree
 *
 *  @return number of times a selected node was from a cut off subtree
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_Longint SCIPgetNDelayedCutoffs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDelayedCutoffs", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->ndelayedcutoffs;
}

/** gets total number of LPs solved so far
 *
 *  @return the total number of LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_Longint SCIPgetNLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPs", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nlps;
}

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm
 *
 *  @return the total number of iterations used so far in primal and dual simplex and barrier algorithm
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nlpiterations;
}

/** gets number of active non-zeros in the current transformed problem
 *
 *  @return the number of active non-zeros in the current transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
 */
SCIP_Longint SCIPgetNNZs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNZs", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nnz;
}

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node
 *
 *  @return the total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNRootLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNRootLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nrootlpiterations;
}

/** gets total number of iterations used in primal and dual simplex and barrier algorithm for the first LP at the root
 *  node
 *
 *  @return the total number of iterations used in primal and dual simplex and barrier algorithm for the first root LP
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNRootFirstLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNRootFirstLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nrootfirstlpiterations;
}

/** gets total number of primal LPs solved so far
 *
 *  @return the total number of primal LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNPrimalLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrimalLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nprimallps;
}

/** gets total number of iterations used so far in primal simplex
 *
 *  @return total number of iterations used so far in primal simplex
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNPrimalLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrimalLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nprimallpiterations;
}

/** gets total number of dual LPs solved so far
 *
 *  @return the total number of dual LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNDualLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDualLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nduallps;
}

/** gets total number of iterations used so far in dual simplex
 *
 *  @return the total number of iterations used so far in dual simplex
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNDualLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDualLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nduallpiterations;
}

/** gets total number of barrier LPs solved so far
 *
 *  @return the total number of barrier LPs solved so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNBarrierLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNBarrierLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nbarrierlps;
}

/** gets total number of iterations used so far in barrier algorithm
 *
 *  @return the total number of iterations used so far in barrier algorithm
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNBarrierLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNBarrierLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nbarrierlpiterations;
}

/** gets total number of LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNResolveLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelps + scip->stat->ndualresolvelps;
}

/** gets total number of simplex iterations used so far in primal and dual simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in primal and dual simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNResolveLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelpiterations + scip->stat->ndualresolvelpiterations;
}

/** gets total number of primal LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of primal LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNPrimalResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrimalResolveLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelps;
}

/** gets total number of simplex iterations used so far in primal simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in primal simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNPrimalResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPrimalResolveLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelpiterations;
}

/** gets total number of dual LPs solved so far that were resolved from an advanced start basis
 *
 *  @return the total number of dual LPs solved so far that were resolved from an advanced start basis
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNDualResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDualResolveLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->ndualresolvelps;
}

/** gets total number of simplex iterations used so far in dual simplex calls where an advanced start basis
 *  was available
 *
 *  @return the total number of simplex iterations used so far in dual simplex calls where an advanced start
 *          basis was available
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNDualResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDualResolveLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->ndualresolvelpiterations;
}

/** gets total number of LPs solved so far for node relaxations
 *
 *  @return the total number of LPs solved so far for node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNNodeLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNodeLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nnodelps;
}

/** gets total number of LPs solved in 0 iterations for node relaxations
 *
 *  @return the total number of LPs solved with 0 iteratins for node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNNodeZeroIterationLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNodeZeroIterationLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nnodezeroitlps;
}

/** gets total number of simplex iterations used so far for node relaxations
 *
 *  @return the total number of simplex iterations used so far for node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNNodeLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNodeLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nnodelpiterations;
}

/** gets total number of LPs solved so far for initial LP in node relaxations
 *
 *  @return the total number of LPs solved so far for initial LP in node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNNodeInitLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNodeInitLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->ninitlps;
}

/** gets total number of simplex iterations used so far for initial LP in node relaxations
 *
 *  @return the total number of simplex iterations used so far for initial LP in node relaxations
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNNodeInitLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNNodeInitLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->ninitlpiterations;
}

/** gets total number of LPs solved so far during diving and probing
 *
 *  @return total number of LPs solved so far during diving and probing
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNDivingLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDivingLPs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->ndivinglps;
}

/** gets total number of simplex iterations used so far during diving and probing
 *
 *  @return the total number of simplex iterations used so far during diving and probing
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNDivingLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNDivingLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->ndivinglpiterations;
}

/** gets total number of times, strong branching was called (each call represents solving two LPs)
 *
 *  @return the total number of times, strong branching was called (each call represents solving two LPs)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNStrongbranchs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nstrongbranchs;
}

/** gets total number of simplex iterations used so far in strong branching
 *
 *  @return the total number of simplex iterations used so far in strong branching
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNStrongbranchLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nsblpiterations;
}

/** gets total number of times, strong branching was called at the root node (each call represents solving two LPs)
 *
 *  @return the total number of times, strong branching was called at the root node (each call represents solving two LPs)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNRootStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNRootStrongbranchs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nrootstrongbranchs;
}

/** gets total number of simplex iterations used so far in strong branching at the root node
 *
 *  @return the total number of simplex iterations used so far in strong branching at the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Longint SCIPgetNRootStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNRootStrongbranchLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nrootsblpiterations;
}

/** gets number of pricing rounds performed so far at the current node
 *
 *  @return the number of pricing rounds performed so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNPriceRounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPriceRounds", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->stat->npricerounds;
}

/** get current number of variables in the pricing store
 *
 *  @return the current number of variables in the pricing store
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNPricevars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPricevars", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->pricestore == NULL ? 0 : SCIPpricestoreGetNVars(scip->pricestore);
}

/** get total number of pricing variables found so far
 *
 *  @return the total number of pricing variables found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNPricevarsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPricevarsFound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->pricestore == NULL ? 0 : SCIPpricestoreGetNVarsFound(scip->pricestore);
}

/** get total number of pricing variables applied to the LPs
 *
 *  @return the total number of pricing variables applied to the LPs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNPricevarsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNPricevarsApplied", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->pricestore == NULL ? 0 : SCIPpricestoreGetNVarsApplied(scip->pricestore);
}

/** gets number of separation rounds performed so far at the current node
 *
 *  @return the number of separation rounds performed so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNSepaRounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNSepaRounds", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->stat->nseparounds;
}

/** get total number of cuts added to the sepastore so far; this includes global cuts from the cut pool as often as they are separated
 *
 *  @return the total number of cuts added to the sepastore so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNCutsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNCutsFound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->sepastore == NULL ? 0 : SCIPsepastoreGetNCutsAdded(scip->sepastore);
}

/** get number of cuts found so far in current separation round
 *
 *  @return the number of cuts found so far in current separation round
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNCutsFoundRound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNCutsFoundRound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->sepastore == NULL ? 0 : SCIPsepastoreGetNCutsFoundRound(scip->sepastore);
}

/** get total number of cuts applied to the LPs
 *
 *  @return the total number of cuts applied to the LPs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNCutsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNCutsApplied", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->sepastore == NULL ? 0 : SCIPsepastoreGetNCutsApplied(scip->sepastore);
}

/** get total number of constraints found in conflict analysis (conflict, reconvergence constraints, and dual proofs)
 *
 *  @return the total number of constraints found in conflict analysis (conflict, reconvergence constraints, and dual proofs)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetNConflictConssFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNConflictConssFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->conflict == NULL ? 0 : (SCIPconflictGetNPropConflictConss(scip->conflict)
      + SCIPconflictGetNPropReconvergenceConss(scip->conflict)
      + SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict)
      + SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict)
      + SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict)
      + SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict)
      + SCIPconflictGetNStrongbranchConflictConss(scip->conflict)
      + SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict)
      + SCIPconflictGetNPseudoConflictConss(scip->conflict)
      + SCIPconflictGetNPseudoReconvergenceConss(scip->conflict)
      + SCIPconflictGetNDualproofsBndGlobal(scip->conflict)
      + SCIPconflictGetNDualproofsInfGlobal(scip->conflict));
}

/** get number of conflict constraints found so far at the current node
 *
 *  @return the number of conflict constraints found so far at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetNConflictConssFoundNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNConflictConssFoundNode", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->conflict == NULL ? 0 : SCIPconflictGetNConflicts(scip->conflict);
}

/** get total number of conflict constraints added to the problem
 *
 *  @return the total number of conflict constraints added to the problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetNConflictConssApplied(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNConflictConssApplied", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->conflict == NULL ? 0 : SCIPconflictGetNAppliedConss(scip->conflict);
}

/** get total number of dual proof constraints added to the problem
 *
 *  @return the total number of dual proof constraints added to the problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetNConflictDualproofsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNConflictDualproofsApplied", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->conflict == NULL ? 0 : (SCIPconflictGetNDualproofsInfSuccess(scip->conflict) +
      SCIPconflictGetNDualproofsBndSuccess(scip->conflict));
}

/** gets maximal depth of all processed nodes in current branch and bound run (excluding probing nodes)
 *
 *  @return the maximal depth of all processed nodes in current branch and bound run (excluding probing nodes)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetMaxDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetMaxDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->maxdepth;
}

/** gets maximal depth of all processed nodes over all branch and bound runs
 *
 *  @return the maximal depth of all processed nodes over all branch and bound runs
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
int SCIPgetMaxTotalDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetMaxTotalDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->maxtotaldepth;
}

/** gets total number of backtracks, i.e. number of times, the new node was selected from the leaves queue
 *
 *  @return the total number of backtracks, i.e. number of times, the new node was selected from the leaves queue
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetNBacktracks(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNBacktracks", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nbacktracks;
}

/** gets total number of active constraints at the current node
 *
 *  @return the total number of active constraints at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNActiveConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNActiveConss", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->stat->nactiveconss;
}

/** gets total number of enabled constraints at the current node
 *
 *  @return the total number of enabled constraints at the current node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
int SCIPgetNEnabledConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNEnabledConss", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   return scip->stat->nenabledconss;
}

/** gets average dual bound of all unprocessed nodes for original problem
 *
 *  @return the average dual bound of all unprocessed nodes for original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set,
	 SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound));
}

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem
 *
 *  @return the average lower (dual) bound of all unprocessed nodes in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound);
}

/** gets global dual bound
 *
 *  @return the global dual bound
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Real SCIPgetDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDualbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* in case we are in presolving we use the stored dual bound if it exits */
   if( scip->set->stage <= SCIP_STAGE_INITSOLVE && scip->transprob->dualbound < SCIP_INVALID )
      return scip->transprob->dualbound;

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPgetLowerbound(scip));
}

/** gets global lower (dual) bound in transformed problem
 *
 *  @return the global lower (dual) bound in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLowerbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->set->stage <= SCIP_STAGE_INITSOLVE )
      return -SCIPinfinity(scip);
   else if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD || SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED )
   {
      /* in case we could not prove whether the problem is unbounded or infeasible, we want to terminate with lower
       * bound = -inf instead of lower bound = upper bound = +inf also in case we prove that the problem is unbounded,
       * it seems to make sense to return with lower bound = -inf, since -infinity is the only valid lower bound
       */
      return -SCIPinfinity(scip);
   }
   else if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE )
   {
      /* SCIPtreeGetLowerbound() should return +inf in the case of infeasibility, but when infeasibility is detected
       * during presolving this does not seem to be the case; hence, we treat this case explicitly
       */
      return SCIPinfinity(scip);
   }
   else
   {
      SCIP_Real treelowerbound;

      /* it may happen that the remaining tree is empty or all open nodes have a lower bound above the cutoff bound, but
       * have not yet been cut off, e.g., when the user calls SCIPgetDualbound() in some event handler; in this case,
       * the global lower bound is given by the upper bound value
       */
      treelowerbound = SCIPtreeGetLowerbound(scip->tree, scip->set);

      if( treelowerbound < scip->primal->upperbound)
         return treelowerbound;
      else
         return scip->primal->upperbound;
   }
}

/** gets dual bound of the root node for the original problem
 *
 *  @return the dual bound of the root node for the original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetDualboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDualboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPgetLowerboundRoot(scip));
}

/** gets lower (dual) bound in transformed problem of the root node
 *
 *  @return the lower (dual) bound in transformed problem of the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetLowerboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetLowerboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->rootlowerbound;
}

/** gets dual bound for the original problem obtained by the first LP solve at the root node
 *
 *  @return the dual bound for the original problem of the first LP solve at the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetFirstLPDualboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetFirstLPDualboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->firstlpdualbound;
}

/** gets lower (dual) bound in transformed problem obtained by the first LP solve at the root node
 *
 *  @return the lower (dual) bound in transformed problem obtained by first LP solve at the root node
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetFirstLPLowerboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetFirstLPLowerboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->stat->firstlpdualbound == SCIP_INVALID ) /*lint !e777*/
      return -SCIPinfinity(scip);
   else
      return SCIPprobInternObjval(scip->transprob, scip->origprob, scip->set, scip->stat->firstlpdualbound);
}

/** gets the primal bound of the very first solution
 *
 * @return the primal bound of the very first solution
 */
SCIP_Real SCIPgetFirstPrimalBound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return scip->stat->firstprimalbound;
}

/** gets global primal bound (objective value of best solution or user objective limit) for the original problem
 *
 *  @return the global primal bound (objective value of best solution or user objective limit) for the original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Real SCIPgetPrimalbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPrimalbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->origprob, scip->set, SCIPgetUpperbound(scip));
}

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit)
 *
 *  @return the global upper (primal) bound in transformed problem (objective value of best solution or user objective limit)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Real SCIPgetUpperbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetUpperbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED )
      return -SCIPinfinity(scip);
   else
      return scip->primal->upperbound;
}

/** gets global cutoff bound in transformed problem: a sub problem with lower bound larger than the cutoff
 *  cannot contain a better feasible solution; usually, this bound is equal to the upper bound, but if the
 *  objective value is always integral, the cutoff bound is (nearly) one less than the upper bound;
 *  additionally, due to objective function domain propagation, the cutoff bound can be further reduced
 *
 *  @return global cutoff bound in transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Real SCIPgetCutoffbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetCutoffbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->primal->cutoffbound;
}

/** updates the cutoff bound
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note using this method in the solving stage can lead to an erroneous SCIP solving status; in particular,
 *        if a solution not respecting the cutoff bound was found before installing a cutoff bound which
 *        renders the remaining problem infeasible, this solution may be reported as optimal
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  @note the given cutoff bound has to better or equal to known one (SCIPgetCutoffbound())
 *  @note a given cutoff bound is also used for updating the objective limit, if possible
 */
SCIP_RETCODE SCIPupdateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             cutoffbound         /**< new cutoff bound */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPupdateCutoffbound", FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   assert(cutoffbound <= SCIPgetCutoffbound(scip));

   SCIP_CALL( SCIPprimalSetCutoffbound(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
         scip->eventfilter, scip->transprob, scip->origprob, scip->tree, scip->reopt, scip->lp, cutoffbound, FALSE) );

   return SCIP_OKAY;
}


/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 *
 *  @return TRUE if the current primal bound is justified with a feasible primal solution, otherwise FALSE
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Bool SCIPisPrimalboundSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPisPrimalboundSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprimalUpperboundIsSol(scip->primal, scip->set, scip->transprob, scip->origprob);
}

/** gets current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 *
 *  @return the current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetGap(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetGap", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* in case we could not prove whether the problem is unbounded or infeasible, we want to terminate with gap = +inf;
    * if the problem was proven to be unbounded or proven to be infeasible we return gap = 0
    */
   if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD )
      return SCIPsetInfinity(scip->set);
   else if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED )
      return 0.0;

   /* the lowerbound is infinity, but SCIP may not have updated the status; in this case, the problem was already solved
    * so we return gap = 0
    */
   if( SCIPsetIsInfinity(scip->set, SCIPgetLowerbound(scip)) )
      return 0.0;

   return SCIPcomputeGap(SCIPsetEpsilon(scip->set), SCIPsetInfinity(scip->set), SCIPgetPrimalbound(scip), SCIPgetDualbound(scip));
}

/** gets current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 *
 *  @return current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetTransGap(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetTransGap", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* in case we could not prove whether the problem is unbounded or infeasible, we want to terminate with gap = +inf;
    * if the problem was proven to be unbounded or proven to be infeasible we return gap = 0
    */
   if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD )
      return SCIPsetInfinity(scip->set);
   else if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED )
      return 0.0;

   /* the lowerbound is infinity, but SCIP may not have updated the status; in this case, the problem was already solved
    * so we return gap = 0
    */
   if( SCIPsetIsInfinity(scip->set, SCIPgetLowerbound(scip)) )
      return 0.0;

   return SCIPcomputeGap(SCIPsetEpsilon(scip->set), SCIPsetInfinity(scip->set), SCIPgetUpperbound(scip), SCIPgetLowerbound(scip));
}

/** gets number of feasible primal solutions found so far
 *
 *  @return the number of feasible primal solutions found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetNSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNSolsFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->primal->nsolsfound;
}

/** gets number of feasible primal solutions respecting the objective limit found so far
 *
 *  @return the number of feasible primal solutions respecting the objective limit found so far
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
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
 */
SCIP_Longint SCIPgetNLimSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED)
      return 0;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNLimSolsFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->primal->nlimsolsfound;
}

/** gets number of feasible primal solutions found so far, that improved the primal bound at the time they were found
 *
 *  @return the number of feasible primal solutions found so far, that improved the primal bound at the time they were found
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 */
SCIP_Longint SCIPgetNBestSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNBestSolsFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->primal->nbestsolsfound;
}

/** gets the average pseudo cost value for the given direction over all variables
 *
 *  @return the average pseudo cost value for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocost(scip->stat->glbhistory, solvaldelta);
}

/** gets the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgPseudocostCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocost(scip->stat->glbhistorycrun, solvaldelta);
}

/** gets the average number of pseudo cost updates for the given direction over all variables
 *
 *  @return the average number of pseudo cost updates for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgPseudocostCount", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocostCount(scip->stat->glbhistory, dir)
      / MAX(scip->transprob->nbinvars + scip->transprob->nintvars, 1);
}

/** gets the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgPseudocostCountCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocostCount(scip->stat->glbhistorycrun, dir)
      / MAX(scip->transprob->nbinvars + scip->transprob->nintvars, 1);
}

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5
 *
 *  @return the average pseudo cost score value over all variables, assuming a fractionality of 0.5
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgPseudocostScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgPseudocostScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   pscostdown = SCIPhistoryGetPseudocost(scip->stat->glbhistory, -0.5);
   pscostup = SCIPhistoryGetPseudocost(scip->stat->glbhistory, +0.5);

   return SCIPbranchGetScore(scip->set, NULL, pscostdown, pscostup);
}

/** gets the average discounted pseudo cost score value over all variables, assuming a fractionality of 0.5
 *
 *  This combines both pscost and ancpscost fields.
 *
 *  @return the average discounted pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  combining both pscost and ancpscost fields
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgDPseudocostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             discountfac         /**< discount factor for discounted pseudocost */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgDPseudocostScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   pscostdown = SCIPhistoryGetPseudocost(scip->stat->glbhistory, -0.5)
               + discountfac * SCIPhistoryGetAncPseudocost(scip->stat->glbhistory, -0.5);
   pscostup = SCIPhistoryGetPseudocost(scip->stat->glbhistory, +0.5)
            + discountfac * SCIPhistoryGetAncPseudocost(scip->stat->glbhistory, +0.5);
   pscostdown /= (1 + discountfac);
   pscostup /= (1 + discountfac);

   return SCIPbranchGetScore(scip->set, NULL, pscostdown, pscostup);
}

/** returns the variance of pseudo costs for all variables in the requested direction
 *
 *  @return the variance of pseudo costs for all variables in the requested direction
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetPseudocostVariance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        branchdir,          /**< the branching direction, up or down */
   SCIP_Bool             onlycurrentrun      /**< use only history of current run? */
   )
{
   SCIP_HISTORY* history;

   assert(scip != NULL);
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPseudocostVariance", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   history = (onlycurrentrun ? scip->stat->glbhistorycrun : scip->stat->glbhistory);
   assert(history != NULL);

   return SCIPhistoryGetPseudocostVariance(history, branchdir);
}

/** gets the number of pseudo cost updates for the given direction over all variables
 *
 *  @return the number of pseudo cost updates for the given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir,                /**< branching direction (downwards, or upwards) */
   SCIP_Bool             onlycurrentrun      /**< use only history of current run? */
   )
{
   SCIP_HISTORY* history;

   assert(scip != NULL);
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetPseudocostCount", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   history = (onlycurrentrun ? scip->stat->glbhistorycrun : scip->stat->glbhistory);

   return SCIPhistoryGetPseudocostCount(history, dir);
}

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 *
 *  @return the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgPseudocostScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgPseudocostScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   pscostdown = SCIPhistoryGetPseudocost(scip->stat->glbhistorycrun, -0.5);
   pscostup = SCIPhistoryGetPseudocost(scip->stat->glbhistorycrun, +0.5);

   return SCIPbranchGetScore(scip->set, NULL, pscostdown, pscostup);
}

/** gets the average conflict score value over all variables
 *
 *  @return the average conflict score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgConflictScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictscoredown;
   SCIP_Real conflictscoreup;
   SCIP_Real scale;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgConflictScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   scale = scip->transprob->nvars * scip->stat->vsidsweight;
   conflictscoredown = SCIPhistoryGetVSIDS(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) / scale;
   conflictscoreup = SCIPhistoryGetVSIDS(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) / scale;

   return SCIPbranchGetScore(scip->set, NULL, conflictscoredown, conflictscoreup);
}

/** gets the average conflict score value over all variables, only using the conflict score information of the current run
 *
 *  @return the average conflict score value over all variables, only using the conflict score information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgConflictScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictscoredown;
   SCIP_Real conflictscoreup;
   SCIP_Real scale;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgConflictScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   scale = scip->transprob->nvars * scip->stat->vsidsweight;
   conflictscoredown = SCIPhistoryGetVSIDS(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS) / scale;
   conflictscoreup = SCIPhistoryGetVSIDS(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS) / scale;

   return SCIPbranchGetScore(scip->set, NULL, conflictscoredown, conflictscoreup);
}

/** gets the average inference score value over all variables
 *
 *  @return the average inference score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgConflictlengthScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictlengthdown;
   SCIP_Real conflictlengthup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgConflictlengthScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   conflictlengthdown = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   conflictlengthup = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, conflictlengthdown, conflictlengthup);
}

/** gets the average conflictlength score value over all variables, only using the conflictlength information of the
 *  current run
 *
 *  @return the average conflictlength score value over all variables, only using the conflictlength information of the
 *          current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgConflictlengthScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictlengthdown;
   SCIP_Real conflictlengthup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgConflictlengthScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   conflictlengthdown = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   conflictlengthup = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, conflictlengthdown, conflictlengthup);
}

/** returns the average number of inferences found after branching in given direction over all variables
 *
 *  @return the average number of inferences found after branching in given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgInferences", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetAvgInferences(scip->stat->glbhistory, dir);
}

/** returns the average number of inferences found after branching in given direction over all variables,
 *  only using the inference information of the current run
 *
 *  @return the average number of inferences found after branching in given direction over all variables,
 *          only using the inference information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgInferencesCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, dir);
}

/** gets the average inference score value over all variables
 *
 *  @return the average inference score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgInferenceScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real inferencesdown;
   SCIP_Real inferencesup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgInferenceScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   inferencesdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   inferencesup = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, inferencesdown, inferencesup);
}

/** gets the average inference score value over all variables, only using the inference information of the
 *  current run
 *
 *  @return the average inference score value over all variables, only using the inference information of the
 *          current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgInferenceScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real inferencesdown;
   SCIP_Real inferencesup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgInferenceScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   inferencesdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   inferencesup = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, inferencesdown, inferencesup);
}

/** returns the average number of cutoffs found after branching in given direction over all variables
 *
 *  @return the average number of cutoffs found after branching in given direction over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgCutoffs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetAvgCutoffs(scip->stat->glbhistory, dir);
}

/** returns the average number of cutoffs found after branching in given direction over all variables,
 *  only using the cutoff information of the current run
 *
 *  @return the average number of cutoffs found after branching in given direction over all variables,
 *          only using the cutoff information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgCutoffsCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetAvgCutoffs(scip->stat->glbhistorycrun, dir);
}

/** gets the average cutoff score value over all variables
 *
 *  @return the average cutoff score value over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgCutoffScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real cutoffsdown;
   SCIP_Real cutoffsup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgCutoffScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   cutoffsdown = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffsup = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, cutoffsdown, cutoffsup);
}

/** gets the average cutoff score value over all variables, only using the cutoff score information of the current run
 *
 *  @return the average cutoff score value over all variables, only using the cutoff score information of the current run
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgCutoffScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real cutoffsdown;
   SCIP_Real cutoffsup;

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   cutoffsdown = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffsup = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, cutoffsdown, cutoffsup);
}

/** increases the average normalized efficacy of a GMI cut over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPincAvgGMIeff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             gmieff              /**< average normalized GMI cut efficacy over all variables */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPincAvgGMIeff", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPhistoryIncGMIeffSum(scip->stat->glbhistory, gmieff);
}

/** returns the average normalized efficacy of a GMI cut over all variables
 *
 *  @return the average normalized efficacy of a GMI cut over all variables
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetAvgGMIeff(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetAvgGMIeff", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPhistoryGetAvgGMIeff(scip->stat->glbhistory);
}

/** computes a deterministic measure of time from statistics
 *
 *  @return the deterministic  time
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_Real SCIPgetDeterministicTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
/* TODO:    SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetDeterministicTime", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) ); */
   if(scip->stat == NULL)
      return 0.0;

   return 1e-6 * scip->stat->nnz * (
          0.00328285264101 * scip->stat->nprimalresolvelpiterations +
          0.00531625104146 * scip->stat->ndualresolvelpiterations +
          0.000738719124051 * scip->stat->nprobboundchgs +
          0.0011123144764 * scip->stat->nisstoppedcalls );
}

/** outputs problem to file stream */
static
SCIP_RETCODE printProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROB*            prob,               /**< problem data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format) */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RESULT result;
   int i;
   assert(scip != NULL);
   assert(prob != NULL);

   /* try all readers until one could read the file */
   result = SCIP_DIDNOTRUN;
   for( i = 0; i < scip->set->nreaders && result == SCIP_DIDNOTRUN; ++i )
   {
      SCIP_RETCODE retcode;

      if( extension != NULL )
         retcode = SCIPreaderWrite(scip->set->readers[i], prob, scip->set, file, extension, genericnames, &result);
      else
         retcode = SCIPreaderWrite(scip->set->readers[i], prob, scip->set, file, "cip", genericnames, &result);

      /* check for reader errors */
      if( retcode == SCIP_WRITEERROR )
         return retcode;

      SCIP_CALL( retcode );
   }

   switch( result )
   {
   case SCIP_DIDNOTRUN:
      return SCIP_PLUGINNOTFOUND;

   case SCIP_SUCCESS:
      return SCIP_OKAY;

   default:
      assert(i < scip->set->nreaders);
      SCIPerrorMessage("invalid result code <%d> from reader <%s> writing <%s> format\n",
         result, SCIPreaderGetName(scip->set->readers[i]), extension);
      return SCIP_READERROR;
   }  /*lint !e788*/
}

/** outputs original problem to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_RETCODE SCIPprintOrigProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintOrigProblem", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert(scip != NULL);
   assert( scip->origprob != NULL );

   retcode = printProblem(scip, scip->origprob, file, extension, genericnames);

   /* check for write errors */
   if( retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** outputs transformed problem of the current node to file stream
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
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
SCIP_RETCODE SCIPprintTransProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintTransProblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   assert(scip != NULL);
   assert(scip->transprob != NULL );

   retcode = printProblem(scip, scip->transprob, file, extension, genericnames);

   /* check for write errors */
   if( retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** outputs status statistics
 *
 *  @note If limits have been changed between the solution and the call to this function, the status is recomputed and
 *        thus may correspond to the original status.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintStatusStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintStatusStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "SCIP Status        : ");
   SCIP_CALL_ABORT( SCIPprintStage(scip, file) );
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");
}

/** collects status statistics in a SCIP_DATATREE object
 *
 *  This function sets: 
 *   - status: the current status of the solver
 *   - info: info about the keys and values stored in the datatree
 *
 *  @note If limits have been changed between the solution and the call to this function, the status is recomputed and
 *        thus may correspond to the original status.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcollectStatusStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectStatusStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPinsertDatatreeString(scip, datatree, "status", SCIPstatusName(SCIPgetStatus(scip))) );

   return SCIP_OKAY;
}

/** outputs statistics for original problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintOrigProblemStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintOrigProblemStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Original Problem   :\n");
   SCIPprobPrintStatistics(scip->origprob, scip->set, scip->messagehdlr, file);
}

/** collects statistics for original problem in a SCIP_DATATREE object
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcollectOrigProblemStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectOrigProblemStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPprobCollectStatistics(scip->origprob, scip->mem->setmem, scip->set, datatree) );

   return SCIP_OKAY;
}

/** outputs statistics for transformed problem
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintTransProblemStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintTransProblemStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Presolved Problem  :\n");
   SCIPprobPrintStatistics(scip->transprob, scip->set, scip->messagehdlr, file);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Nonzeros         : %" SCIP_LONGINT_FORMAT " constraint, %" SCIP_LONGINT_FORMAT " clique table\n",
         scip->stat->nnz, SCIPcliquetableGetNEntries(scip->cliquetable));
}

/** collects statistics for transformed problem in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectTransProblemStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectTransProblemStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Collect problem statistics */
   SCIP_CALL( SCIPprobCollectStatistics(scip->transprob, scip->mem->setmem, scip->set, datatree) );

   /* Collect additional statistics */
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "constraint_nonzeros", scip->stat->nnz) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "clique_table_nonzeros", SCIPcliquetableGetNEntries(scip->cliquetable)) );

   return SCIP_OKAY;
}

/** outputs presolver statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintPresolverStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintPresolverStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Presolvers         :   ExecTime  SetupTime  Calls  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs\n");

   /* sort presolvers w.r.t. their name */
   SCIPsetSortPresolsName(scip->set);

   /* presolver statistics */
   for( i = 0; i < scip->set->npresols; ++i )
   {
      SCIP_PRESOL* presol;
      presol = scip->set->presols[i];
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s:", SCIPpresolGetName(presol));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10.2f %6d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
         SCIPpresolGetTime(presol),
         SCIPpresolGetSetupTime(presol),
         SCIPpresolGetNCalls(presol),
         SCIPpresolGetNFixedVars(presol),
         SCIPpresolGetNAggrVars(presol),
         SCIPpresolGetNChgVarTypes(presol),
         SCIPpresolGetNChgBds(presol),
         SCIPpresolGetNAddHoles(presol),
         SCIPpresolGetNDelConss(presol),
         SCIPpresolGetNAddConss(presol),
         SCIPpresolGetNChgSides(presol),
         SCIPpresolGetNChgCoefs(presol));
   }

   /* sort propagators w.r.t. their name */
   SCIPsetSortPropsName(scip->set);

   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop;
      prop = scip->set->props[i];
      if( SCIPpropDoesPresolve(prop) )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s:", SCIPpropGetName(prop));
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10.2f %6d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
            SCIPpropGetPresolTime(prop),
	    SCIPpropGetSetupTime(prop),
            SCIPpropGetNPresolCalls(prop),
            SCIPpropGetNFixedVars(prop),
            SCIPpropGetNAggrVars(prop),
            SCIPpropGetNChgVarTypes(prop),
            SCIPpropGetNChgBds(prop),
            SCIPpropGetNAddHoles(prop),
            SCIPpropGetNDelConss(prop),
            SCIPpropGetNAddConss(prop),
            SCIPpropGetNChgSides(prop),
            SCIPpropGetNChgCoefs(prop));
      }
   }

   /* constraint handler presolving methods statistics */
   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      int maxnactiveconss;

      conshdlr = scip->set->conshdlrs[i];
      maxnactiveconss = SCIPconshdlrGetMaxNActiveConss(conshdlr);
      if( SCIPconshdlrDoesPresolve(conshdlr)
         && (maxnactiveconss > 0 || !SCIPconshdlrNeedsCons(conshdlr)
            || SCIPconshdlrGetNFixedVars(conshdlr) > 0
            || SCIPconshdlrGetNAggrVars(conshdlr) > 0
            || SCIPconshdlrGetNChgVarTypes(conshdlr) > 0
            || SCIPconshdlrGetNChgBds(conshdlr) > 0
            || SCIPconshdlrGetNAddHoles(conshdlr) > 0
            || SCIPconshdlrGetNDelConss(conshdlr) > 0
            || SCIPconshdlrGetNAddConss(conshdlr) > 0
            || SCIPconshdlrGetNChgSides(conshdlr) > 0
            || SCIPconshdlrGetNChgCoefs(conshdlr) > 0
            || SCIPconshdlrGetNUpgdConss(conshdlr) > 0) )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10.2f %6d %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
            SCIPconshdlrGetPresolTime(conshdlr),
            SCIPconshdlrGetSetupTime(conshdlr),
            SCIPconshdlrGetNPresolCalls(conshdlr),
            SCIPconshdlrGetNFixedVars(conshdlr),
            SCIPconshdlrGetNAggrVars(conshdlr),
            SCIPconshdlrGetNChgVarTypes(conshdlr),
            SCIPconshdlrGetNChgBds(conshdlr),
            SCIPconshdlrGetNAddHoles(conshdlr),
            SCIPconshdlrGetNDelConss(conshdlr),
            SCIPconshdlrGetNAddConss(conshdlr),
            SCIPconshdlrGetNChgSides(conshdlr),
            SCIPconshdlrGetNChgCoefs(conshdlr));
      }
   }

   /* root node bound changes */
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  root node        :          -          -      - %10d          -          - %10d          -          -          -          -          -\n",
      scip->stat->nrootintfixings, scip->stat->nrootboundchgs);
}

/** collects presolver statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectPresolverStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* plugins;
   SCIP_DATATREE* rootdata;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectPresolverStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Sort presolvers by name */
   SCIPsetSortPresolsName(scip->set);

   SCIP_CALL( SCIPcreateDatatreeInTree( scip, datatree, &plugins, "plugins", scip->set->npresols + scip->set->nprops + scip->set->nconshdlrs) );

   /* Collect presolver statistics */
   for( i = 0; i < scip->set->npresols; ++i )
   {
      SCIP_PRESOL* presol = scip->set->presols[i];
      SCIP_DATATREE* presoldata;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, plugins, &presoldata, SCIPpresolGetName(presol), 13) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, presoldata, "description", SCIPpresolGetDesc(presol)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, presoldata, "exec_time", SCIPpresolGetTime(presol)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, presoldata, "setup_time", SCIPpresolGetSetupTime(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "calls", SCIPpresolGetNCalls(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "fixed_vars", SCIPpresolGetNFixedVars(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "aggregated_vars", SCIPpresolGetNAggrVars(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "changed_var_types", SCIPpresolGetNChgVarTypes(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "changed_bounds", SCIPpresolGetNChgBds(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "added_holes", SCIPpresolGetNAddHoles(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "deleted_constraints", SCIPpresolGetNDelConss(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "added_constraints", SCIPpresolGetNAddConss(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "changed_sides", SCIPpresolGetNChgSides(presol)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, presoldata, "changed_coefficients", SCIPpresolGetNChgCoefs(presol)) );
   }

   /* Sort propagators by name */
   SCIPsetSortPropsName(scip->set);

   /* Collect propagator statistics */
   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop = scip->set->props[i];
      SCIP_DATATREE* propdata;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, plugins, &propdata, SCIPpropGetName(prop), 15) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, propdata, "description", SCIPpropGetDesc(prop)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "exec_time", SCIPpropGetTime(prop)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "setup_time", SCIPpropGetSetupTime(prop)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, propdata, "calls", SCIPpropGetNCalls(prop)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, propdata, "cutoffs", SCIPpropGetNCutoffs(prop)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, propdata, "domain_reductions", SCIPpropGetNDomredsFound(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "fixed_vars", SCIPpropGetNFixedVars(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "aggregated_vars", SCIPpropGetNAggrVars(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "changed_var_types", SCIPpropGetNChgVarTypes(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "changed_bounds", SCIPpropGetNChgBds(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "added_holes", SCIPpropGetNAddHoles(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "deleted_constraints", SCIPpropGetNDelConss(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "added_constraints", SCIPpropGetNAddConss(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "changed_sides", SCIPpropGetNChgSides(prop)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, propdata, "changed_coefficients", SCIPpropGetNChgCoefs(prop)) );
   }

   /* Collect constraint handler presolving methods statistics */
   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr = scip->set->conshdlrs[i];
      if( SCIPconshdlrDoesPresolve(conshdlr) )
      {
         SCIP_DATATREE* conshdlrdata;

         SCIP_CALL( SCIPcreateDatatreeInTree( scip, plugins, &conshdlrdata, SCIPconshdlrGetName(conshdlr), 13) );

         SCIP_CALL( SCIPinsertDatatreeString(scip, conshdlrdata, "description", SCIPconshdlrGetDesc(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrdata, "presol_time", SCIPconshdlrGetPresolTime(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrdata, "setup_time", SCIPconshdlrGetSetupTime(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "presol_calls", SCIPconshdlrGetNPresolCalls(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "fixed_vars", SCIPconshdlrGetNFixedVars(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "aggregated_vars", SCIPconshdlrGetNAggrVars(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "changed_var_types", SCIPconshdlrGetNChgVarTypes(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "changed_bounds", SCIPconshdlrGetNChgBds(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "added_holes", SCIPconshdlrGetNAddHoles(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "deleted_constraints", SCIPconshdlrGetNDelConss(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "added_constraints", SCIPconshdlrGetNAddConss(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "changed_sides", SCIPconshdlrGetNChgSides(conshdlr)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "changed_coefficients", SCIPconshdlrGetNChgCoefs(conshdlr)) );
      }
   }

   /* Collect root node fixings statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &rootdata, "root", 2) );

   SCIP_CALL( SCIPinsertDatatreeInt(scip, rootdata, "int_fixings", scip->stat->nrootintfixings ) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, rootdata, "bound_changes", scip->stat->nrootboundchgs ) );

   return SCIP_OKAY;
}

/** outputs constraint statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintConstraintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintConstraintStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Add maximal number of constraints of the same type? So far this information is not added because of lack of space. */
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoRelax  #EnfoPS    #Check   #ResProp    Cutoffs    DomReds       Cuts    Applied      Conss   Children\n");

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      int startnactiveconss;
      int maxnactiveconss;

      conshdlr = scip->set->conshdlrs[i];
      startnactiveconss = SCIPconshdlrGetStartNActiveConss(conshdlr);
      maxnactiveconss = SCIPconshdlrGetMaxNActiveConss(conshdlr);
      if( maxnactiveconss > 0 || !SCIPconshdlrNeedsCons(conshdlr) )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10d%c%10d %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
            startnactiveconss,
            maxnactiveconss > startnactiveconss ? '+' : ' ',
            maxnactiveconss,
            SCIPconshdlrGetNSepaCalls(conshdlr),
            SCIPconshdlrGetNPropCalls(conshdlr),
            SCIPconshdlrGetNEnfoLPCalls(conshdlr),
            SCIPconshdlrGetNEnfoRelaxCalls(conshdlr),
            SCIPconshdlrGetNEnfoPSCalls(conshdlr),
            SCIPconshdlrGetNCheckCalls(conshdlr),
            SCIPconshdlrGetNRespropCalls(conshdlr),
            SCIPconshdlrGetNCutoffs(conshdlr),
            SCIPconshdlrGetNDomredsFound(conshdlr),
            SCIPconshdlrGetNCutsFound(conshdlr),
            SCIPconshdlrGetNCutsApplied(conshdlr),
            SCIPconshdlrGetNConssFound(conshdlr),
            SCIPconshdlrGetNChildren(conshdlr));
      }
   }
}

/** collects constraint statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectConstraintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* constraints;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectConstraintStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Create a subtree for constraints */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &constraints, "plugins", scip->set->nconshdlrs) );

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr = scip->set->conshdlrs[i];
      SCIP_DATATREE* conshdlrdata;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, constraints, &conshdlrdata, SCIPconshdlrGetName(conshdlr), 16) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, conshdlrdata, "description", SCIPconshdlrGetDesc(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "start_active_constraints", SCIPconshdlrGetStartNActiveConss(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, conshdlrdata, "max_active_constraints", SCIPconshdlrGetMaxNActiveConss(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "separation_calls", SCIPconshdlrGetNSepaCalls(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "propagation_calls", SCIPconshdlrGetNPropCalls(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "enforcement_lp_calls", SCIPconshdlrGetNEnfoLPCalls(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "enforcement_relax_calls", SCIPconshdlrGetNEnfoRelaxCalls(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "enforcement_ps_calls", SCIPconshdlrGetNEnfoPSCalls(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "check_calls", SCIPconshdlrGetNCheckCalls(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "response_propagation_calls", SCIPconshdlrGetNRespropCalls(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "cutoffs", SCIPconshdlrGetNCutoffs(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "domain_reductions", SCIPconshdlrGetNDomredsFound(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "cuts_found", SCIPconshdlrGetNCutsFound(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "cuts_applied", SCIPconshdlrGetNCutsApplied(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "constraints_found", SCIPconshdlrGetNConssFound(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, conshdlrdata, "children_created", SCIPconshdlrGetNChildren(conshdlr)) );
   }

   return SCIP_OKAY;
}

/** outputs constraint timing statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintConstraintTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintConstraintTimingStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Constraint Timings :  TotalTime  SetupTime   Separate  Propagate     EnfoLP     EnfoPS     EnfoRelax   Check    ResProp    SB-Prop\n");

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      int maxnactiveconss;

      conshdlr = scip->set->conshdlrs[i];
      maxnactiveconss = SCIPconshdlrGetMaxNActiveConss(conshdlr);
      if( maxnactiveconss > 0 || !SCIPconshdlrNeedsCons(conshdlr) )
      {
         SCIP_Real totaltime;

         totaltime = SCIPconshdlrGetSepaTime(conshdlr) + SCIPconshdlrGetPropTime(conshdlr)
            + SCIPconshdlrGetStrongBranchPropTime(conshdlr)
            + SCIPconshdlrGetEnfoLPTime(conshdlr)
            + SCIPconshdlrGetEnfoPSTime(conshdlr)
            + SCIPconshdlrGetEnfoRelaxTime(conshdlr)
            + SCIPconshdlrGetCheckTime(conshdlr)
            + SCIPconshdlrGetRespropTime(conshdlr)
	    + SCIPconshdlrGetSetupTime(conshdlr);

         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
            totaltime,
	    SCIPconshdlrGetSetupTime(conshdlr),
            SCIPconshdlrGetSepaTime(conshdlr),
            SCIPconshdlrGetPropTime(conshdlr),
            SCIPconshdlrGetEnfoLPTime(conshdlr),
            SCIPconshdlrGetEnfoPSTime(conshdlr),
            SCIPconshdlrGetEnfoRelaxTime(conshdlr),
            SCIPconshdlrGetCheckTime(conshdlr),
            SCIPconshdlrGetRespropTime(conshdlr),
            SCIPconshdlrGetStrongBranchPropTime(conshdlr));
      }
   }
}

/** collects constraint timing statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectConstraintTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* constrainttimings;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectConstraintTimingStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Create a subtree for constraint timings */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &constrainttimings, "plugins", scip->set->nconshdlrs) );

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr = scip->set->conshdlrs[i];
      SCIP_DATATREE* conshdlrtimingdata;
      SCIP_Real totaltime;

      SCIP_CALL( SCIPcreateDatatreeInTree( scip, constrainttimings, &conshdlrtimingdata, SCIPconshdlrGetName(conshdlr), 10) );

      totaltime = SCIPconshdlrGetSepaTime(conshdlr) + SCIPconshdlrGetPropTime(conshdlr)
         + SCIPconshdlrGetStrongBranchPropTime(conshdlr) + SCIPconshdlrGetEnfoLPTime(conshdlr)
         + SCIPconshdlrGetEnfoPSTime(conshdlr) + SCIPconshdlrGetEnfoRelaxTime(conshdlr)
         + SCIPconshdlrGetCheckTime(conshdlr) + SCIPconshdlrGetRespropTime(conshdlr)
         + SCIPconshdlrGetSetupTime(conshdlr);

      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "total_time", totaltime) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "setup_time", SCIPconshdlrGetSetupTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "separation_time", SCIPconshdlrGetSepaTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "propagation_time", SCIPconshdlrGetPropTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "enforcement_lp_time", SCIPconshdlrGetEnfoLPTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "enforcement_ps_time", SCIPconshdlrGetEnfoPSTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "enforcement_relax_time", SCIPconshdlrGetEnfoRelaxTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "check_time", SCIPconshdlrGetCheckTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "response_propagation_time", SCIPconshdlrGetRespropTime(conshdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conshdlrtimingdata, "strong_branch_propagation_time", SCIPconshdlrGetStrongBranchPropTime(conshdlr)) );
   }

   return SCIP_OKAY;
}

/** outputs propagator statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintPropagatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintPropagatorStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Propagators        : #Propagate   #ResProp    Cutoffs    DomReds\n");

   /* sort propagaters w.r.t. their name */
   SCIPsetSortPropsName(scip->set);

   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop;
      prop = scip->set->props[i];

      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
         SCIPpropGetName(prop),
         SCIPpropGetNCalls(prop),
         SCIPpropGetNRespropCalls(prop),
         SCIPpropGetNCutoffs(prop),
         SCIPpropGetNDomredsFound(prop));
   }

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Propagator Timings :  TotalTime  SetupTime   Presolve  Propagate    ResProp    SB-Prop\n");

   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop;
      SCIP_Real totaltime;

      prop = scip->set->props[i];
      totaltime = SCIPpropGetPresolTime(prop) + SCIPpropGetTime(prop) + SCIPpropGetRespropTime(prop)
         + SCIPpropGetStrongBranchPropTime(prop) + SCIPpropGetSetupTime(prop);

      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s:", SCIPpropGetName(prop));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
         totaltime,
	 SCIPpropGetSetupTime(prop),
	 SCIPpropGetPresolTime(prop),
	 SCIPpropGetTime(prop),
	 SCIPpropGetRespropTime(prop),
	 SCIPpropGetStrongBranchPropTime(prop));
   }
}

/** collects propagator statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectPropagatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* propagators;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectPropagatorStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Create a subtree for propagators */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &propagators, "plugins", scip->set->nprops) );

   /* Collect propagator statistics */
   SCIPsetSortPropsName(scip->set);

   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop = scip->set->props[i];
      SCIP_DATATREE* propdata;
      SCIP_Real totaltime;

      SCIP_CALL( SCIPcreateDatatreeInTree( scip, propagators, &propdata, SCIPpropGetName(prop), 11) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, propdata, "description", SCIPpropGetDesc(prop)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, propdata, "calls", SCIPpropGetNCalls(prop)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, propdata, "resprop_calls", SCIPpropGetNRespropCalls(prop)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, propdata, "cutoffs", SCIPpropGetNCutoffs(prop)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, propdata, "domreds_found", SCIPpropGetNDomredsFound(prop)) );

      totaltime = SCIPpropGetPresolTime(prop) + SCIPpropGetTime(prop) + SCIPpropGetRespropTime(prop)
         + SCIPpropGetStrongBranchPropTime(prop) + SCIPpropGetSetupTime(prop);

      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "total_time", totaltime) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "setup_time", SCIPpropGetSetupTime(prop)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "presolve_time", SCIPpropGetPresolTime(prop)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "propagation_time", SCIPpropGetTime(prop)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "response_propagation_time", SCIPpropGetRespropTime(prop)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, propdata, "strong_branch_propagation_time", SCIPpropGetStrongBranchPropTime(prop)) );
   }

   return SCIP_OKAY;
}

/** outputs conflict statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintConflictStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   char initstoresize[SCIP_MAXSTRLEN];
   char maxstoresize[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintConflictStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->set->conf_maxstoresize == 0 )
   {
      (void)SCIPsnprintf(initstoresize, SCIP_MAXSTRLEN, "inf");
      (void)SCIPsnprintf(maxstoresize, SCIP_MAXSTRLEN, "inf");
   }
   else
   {
      int initsize = SCIPconflictstoreGetInitPoolSize(scip->conflictstore);
      int maxsize = SCIPconflictstoreGetMaxPoolSize(scip->conflictstore);

      if( maxsize == -1 )
      {
         (void)SCIPsnprintf(initstoresize, SCIP_MAXSTRLEN, "--");
         (void)SCIPsnprintf(maxstoresize, SCIP_MAXSTRLEN, "--");
      }
      else
      {
         assert(initsize >= 0);
         assert(maxsize >= 0);

         (void)SCIPsnprintf(initstoresize, SCIP_MAXSTRLEN, "%d", initsize);
         (void)SCIPsnprintf(maxstoresize, SCIP_MAXSTRLEN, "%d", maxsize);
      }
   }
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Conflict Analysis  :       Time      Calls    Success    DomReds  Conflicts   Literals    Reconvs ReconvLits   Dualrays   Nonzeros   LP Iters   (pool size: [%s,%s])\n", initstoresize, maxstoresize);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  propagation      : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          - %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT " %10.1f          -          -          -\n",
      SCIPconflictGetPropTime(scip->conflict),
      SCIPconflictGetNPropCalls(scip->conflict),
      SCIPconflictGetNPropSuccess(scip->conflict),
      SCIPconflictGetNPropConflictConss(scip->conflict),
      SCIPconflictGetNPropConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPropConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPropConflictConss(scip->conflict) : 0,
      SCIPconflictGetNPropReconvergenceConss(scip->conflict),
      SCIPconflictGetNPropReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPropReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPropReconvergenceConss(scip->conflict) : 0);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  infeasible LP    : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          - %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT "\n",
      SCIPconflictGetInfeasibleLPTime(scip->conflict),
      SCIPconflictGetNInfeasibleLPCalls(scip->conflict),
      SCIPconflictGetNInfeasibleLPSuccess(scip->conflict),
      SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict),
      SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNInfeasibleLPConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict) : 0,
      SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict),
      SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNInfeasibleLPReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict) : 0,
      SCIPconflictGetNDualproofsInfSuccess(scip->conflict),
      SCIPconflictGetNDualproofsInfSuccess(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNDualproofsInfNonzeros(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNDualproofsInfSuccess(scip->conflict) : 0,
      SCIPconflictGetNInfeasibleLPIterations(scip->conflict));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  bound exceed. LP : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          - %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT "\n",
      SCIPconflictGetBoundexceedingLPTime(scip->conflict),
      SCIPconflictGetNBoundexceedingLPCalls(scip->conflict),
      SCIPconflictGetNBoundexceedingLPSuccess(scip->conflict),
      SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict),
      SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNBoundexceedingLPConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict) : 0,
      SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict),
      SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNBoundexceedingLPReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict) : 0,
      SCIPconflictGetNDualproofsBndSuccess(scip->conflict),
      SCIPconflictGetNDualproofsBndSuccess(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNDualproofsBndNonzeros(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNDualproofsBndSuccess(scip->conflict) : 0,
      SCIPconflictGetNBoundexceedingLPIterations(scip->conflict));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  strong branching : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          - %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT " %10.1f          -          - %10" SCIP_LONGINT_FORMAT "\n",
      SCIPconflictGetStrongbranchTime(scip->conflict),
      SCIPconflictGetNStrongbranchCalls(scip->conflict),
      SCIPconflictGetNStrongbranchSuccess(scip->conflict),
      SCIPconflictGetNStrongbranchConflictConss(scip->conflict),
      SCIPconflictGetNStrongbranchConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNStrongbranchConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNStrongbranchConflictConss(scip->conflict) : 0,
      SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict),
      SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNStrongbranchReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict) : 0,
      SCIPconflictGetNStrongbranchIterations(scip->conflict));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  pseudo solution  : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          - %10" SCIP_LONGINT_FORMAT " %10.1f %10" SCIP_LONGINT_FORMAT " %10.1f          -          -          -\n",
      SCIPconflictGetPseudoTime(scip->conflict),
      SCIPconflictGetNPseudoCalls(scip->conflict),
      SCIPconflictGetNPseudoSuccess(scip->conflict),
      SCIPconflictGetNPseudoConflictConss(scip->conflict),
      SCIPconflictGetNPseudoConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPseudoConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPseudoConflictConss(scip->conflict) : 0,
      SCIPconflictGetNPseudoReconvergenceConss(scip->conflict),
      SCIPconflictGetNPseudoReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPseudoReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPseudoReconvergenceConss(scip->conflict) : 0);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  applied globally : %10.2f          -          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.1f          -          - %10" SCIP_LONGINT_FORMAT "          -          -\n",
      SCIPconflictGetGlobalApplTime(scip->conflict),
      SCIPconflictGetNGlobalChgBds(scip->conflict),
      SCIPconflictGetNAppliedGlobalConss(scip->conflict),
      SCIPconflictGetNAppliedGlobalConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNAppliedGlobalLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNAppliedGlobalConss(scip->conflict) : 0,
      SCIPconflictGetNDualproofsInfGlobal(scip->conflict) + SCIPconflictGetNDualproofsBndGlobal(scip->conflict));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  applied locally  :          -          -          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.1f          -          - %10" SCIP_LONGINT_FORMAT "          -          -\n",
      SCIPconflictGetNLocalChgBds(scip->conflict),
      SCIPconflictGetNAppliedLocalConss(scip->conflict),
      SCIPconflictGetNAppliedLocalConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNAppliedLocalLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNAppliedLocalConss(scip->conflict) : 0,
      SCIPconflictGetNDualproofsInfLocal(scip->conflict) + SCIPconflictGetNDualproofsBndLocal(scip->conflict));
}

/** collects conflict statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectConflictStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   int initstoresize;
   int maxstoresize;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectConflictStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Collect pool size information */
   if( scip->set->conf_maxstoresize == 0 )
   {
      initstoresize = INT_MAX;
      maxstoresize = INT_MAX;
   }
   else
   {
      initstoresize = SCIPconflictstoreGetInitPoolSize(scip->conflictstore);
      maxstoresize = SCIPconflictstoreGetMaxPoolSize(scip->conflictstore);
   }

   /* if maxstoresize is -1, then we have no sizes to report */
   if( maxstoresize >= 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "pool_size", initstoresize) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "max_pool_size", maxstoresize) );
   }

   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "propagation_time", SCIPconflictGetPropTime(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "propagation_calls", SCIPconflictGetNPropCalls(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "propagation_success", SCIPconflictGetNPropSuccess(scip->conflict)) );

   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "infeasible_lp_time", SCIPconflictGetInfeasibleLPTime(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "infeasible_lp_calls", SCIPconflictGetNInfeasibleLPCalls(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "infeasible_lp_success", SCIPconflictGetNInfeasibleLPSuccess(scip->conflict)) );

   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "bound_exceeding_lp_time", SCIPconflictGetBoundexceedingLPTime(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "bound_exceeding_lp_calls", SCIPconflictGetNBoundexceedingLPCalls(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "bound_exceeding_lp_success", SCIPconflictGetNBoundexceedingLPSuccess(scip->conflict)) );

   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "strong_branching_time", SCIPconflictGetStrongbranchTime(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "strong_branching_calls", SCIPconflictGetNStrongbranchCalls(scip->conflict)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "strong_branching_success", SCIPconflictGetNStrongbranchSuccess(scip->conflict)) );

   return SCIP_OKAY;
}

/** outputs separator statistics
 *
 *  Columns:
 *    - RootCalls: The number of calls that happened at the root.
 *    - FoundCuts: The total number of cuts generated by the separators.
 *      Note: Cutpool-FoundCuts \f$= \sum_{i=1}^nsepas ( Foundcuts_i - DirectAdd_i )\f$.
 *    - ViaPoolAdd: The total number of cuts added to the sepastore from the cutpool.
 *    - DirectAdd: The total number of cuts added directly to the sepastore from the separator.
 *    - Applied: The sum of all cuts from the separator that were applied to the LP.
 *    - ViaPoolApp: The number of cuts that entered the sepastore from the cutpool that were applied to the LP.
 *    - DirectApp: The number of cuts that entered the sepastore directly and were applied to the LP.
 *
 *  The number of cuts ViaPoolAdd + Directly should be equal to the number of cuts Filtered + Forced + Selected in the
 *  cutselector statistics.
 *
 *  @note The following edge case may lead to over or undercounting of statistics: When SCIPapplyCutsProbing() is
 *        called, cuts are counted for the cut selection statistics, but not for the separator statistics.  This
 *        happens, e.g., in the default plugin prop_obbt.c.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintSeparatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintSeparatorStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Separators         :   ExecTime  SetupTime      Calls  RootCalls    Cutoffs    DomReds  FoundCuts ViaPoolAdd  DirectAdd    Applied ViaPoolApp  DirectApp      Conss\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  cut pool         : %10.2f          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          -          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          -          -          -          -          -    (maximal pool size: %10" SCIP_LONGINT_FORMAT")\n",
      SCIPcutpoolGetTime(scip->cutpool),
      SCIPcutpoolGetNCalls(scip->cutpool),
      SCIPcutpoolGetNRootCalls(scip->cutpool),
      SCIPcutpoolGetNCutsFound(scip->cutpool),
      SCIPcutpoolGetNCutsAdded(scip->cutpool),
      SCIPcutpoolGetMaxNCuts(scip->cutpool));

   /* sort separators w.r.t. their name */
   SCIPsetSortSepasName(scip->set);

   for( i = 0; i < scip->set->nsepas; ++i )
   {
      SCIP_SEPA* sepa;

      sepa = scip->set->sepas[i];

      /* only output data for separators without parent separator */
      if( SCIPsepaGetParentsepa(sepa) == NULL )
      {
         /* output data */
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
            SCIPsepaGetName(sepa),
            SCIPsepaGetTime(sepa),
            SCIPsepaGetSetupTime(sepa),
            SCIPsepaGetNCalls(sepa),
            SCIPsepaGetNRootCalls(sepa),
            SCIPsepaGetNCutoffs(sepa),
            SCIPsepaGetNDomredsFound(sepa),
            SCIPsepaGetNCutsFound(sepa),
            SCIPsepaGetNCutsAddedViaPool(sepa),
            SCIPsepaGetNCutsAddedDirect(sepa),
            SCIPsepaGetNCutsApplied(sepa),
            SCIPsepaGetNCutsAppliedViaPool(sepa),
            SCIPsepaGetNCutsAppliedDirect(sepa),
            SCIPsepaGetNConssFound(sepa));

         /* for parent separators search for dependent separators */
         if( SCIPsepaIsParentsepa(sepa) )
         {
            SCIP_SEPA* parentsepa;
            int k;

            for( k = 0; k < scip->set->nsepas; ++k )
            {
               if( k == i )
                  continue;

               parentsepa = SCIPsepaGetParentsepa(scip->set->sepas[k]);
               if( parentsepa != sepa )
                  continue;

               SCIPmessageFPrintInfo(scip->messagehdlr, file, "  > %-15.17s: %10s %10s %10s %10s %10s %10s %10s %10" SCIP_LONGINT_FORMAT" %10" SCIP_LONGINT_FORMAT" %10" SCIP_LONGINT_FORMAT" %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10s\n",
                  SCIPsepaGetName(scip->set->sepas[k]), "-", "-", "-", "-", "-", "-", "-",
                  SCIPsepaGetNCutsAddedViaPool(scip->set->sepas[k]),
                  SCIPsepaGetNCutsAddedDirect(scip->set->sepas[k]),
                  SCIPsepaGetNCutsApplied(scip->set->sepas[k]),
                  SCIPsepaGetNCutsAppliedViaPool(scip->set->sepas[k]),
                  SCIPsepaGetNCutsAppliedDirect(scip->set->sepas[k]), "-");
            }
         }
      }
   }
}

/** collects separator statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectSeparatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* separators;
   SCIP_DATATREE* cutpool;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectSeparatorStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Create a subtree for separators */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &separators, "plugins", scip->set->nsepas + 1) );

   /* Collect statistics for the cut pool */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, separators, &cutpool, "cut_pool", 7) );

   SCIP_CALL( SCIPinsertDatatreeReal(scip, cutpool, "exec_time", SCIPcutpoolGetTime(scip->cutpool)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, cutpool, "calls", SCIPcutpoolGetNCalls(scip->cutpool)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, cutpool, "root_calls", SCIPcutpoolGetNRootCalls(scip->cutpool)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, cutpool, "cuts_found", SCIPcutpoolGetNCutsFound(scip->cutpool)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, cutpool, "cuts_added", SCIPcutpoolGetNCutsAdded(scip->cutpool)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, cutpool, "max_cuts", SCIPcutpoolGetMaxNCuts(scip->cutpool)) );

   /* Sort separators by name */
   SCIPsetSortSepasName(scip->set);

   /* Collect statistics for each separator */
   for( i = 0; i < scip->set->nsepas; ++i )
   {
      SCIP_SEPA* sepa = scip->set->sepas[i];

      /* Only collect data for separators without parent separators */
      if( SCIPsepaGetParentsepa(sepa) == NULL )
      {
         SCIP_DATATREE* sepadata;
         SCIP_CALL( SCIPcreateDatatreeInTree( scip, separators, &sepadata, SCIPsepaGetName(sepa), 14) );

         SCIP_CALL( SCIPinsertDatatreeString(scip, sepadata, "description", SCIPsepaGetDesc(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, sepadata, "exec_time", SCIPsepaGetTime(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, sepadata, "setup_time", SCIPsepaGetSetupTime(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "calls", SCIPsepaGetNCalls(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "root_calls", SCIPsepaGetNRootCalls(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "cutoffs", SCIPsepaGetNCutoffs(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "domreds_found", SCIPsepaGetNDomredsFound(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "cuts_found", SCIPsepaGetNCutsFound(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "cuts_added_via_pool", SCIPsepaGetNCutsAddedViaPool(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "cuts_added_direct", SCIPsepaGetNCutsAddedDirect(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "cuts_applied", SCIPsepaGetNCutsApplied(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "cuts_applied_via_pool", SCIPsepaGetNCutsAppliedViaPool(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "cuts_applied_direct", SCIPsepaGetNCutsAppliedDirect(sepa)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, sepadata, "constraints_found", SCIPsepaGetNConssFound(sepa)) );

         /* for parent separators search for dependent separators */
         if( SCIPsepaIsParentsepa(sepa) )
         {
            SCIP_SEPA* parentsepa;
            int k;

            for( k = 0; k < scip->set->nsepas; ++k )
            {
               if( k == i )
                  continue;

               parentsepa = SCIPsepaGetParentsepa(scip->set->sepas[k]);
               if( parentsepa != sepa )
                  continue;

               SCIP_DATATREE* subsepadata;
               SCIP_CALL( SCIPcreateDatatreeInTree(scip, sepadata, &subsepadata, SCIPsepaGetName(scip->set->sepas[k]), 14) );

               SCIP_CALL( SCIPinsertDatatreeString(scip, subsepadata, "description", SCIPsepaGetDesc(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, subsepadata, "exec_time", SCIPsepaGetTime(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, subsepadata, "setup_time", SCIPsepaGetSetupTime(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "calls", SCIPsepaGetNCalls(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "root_calls", SCIPsepaGetNRootCalls(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "cutoffs", SCIPsepaGetNCutoffs(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "domreds_found", SCIPsepaGetNDomredsFound(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "cuts_found", SCIPsepaGetNCutsFound(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "cuts_added_via_pool", SCIPsepaGetNCutsAddedViaPool(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "cuts_added_direct", SCIPsepaGetNCutsAddedDirect(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "cuts_applied", SCIPsepaGetNCutsApplied(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "cuts_applied_via_pool", SCIPsepaGetNCutsAppliedViaPool(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "cuts_applied_direct", SCIPsepaGetNCutsAppliedDirect(scip->set->sepas[k])) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, subsepadata, "constraints_found", SCIPsepaGetNConssFound(scip->set->sepas[k])) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** outputs cutselector statistics
 *
 *       Filtered = ViaPoolAdd(Separators) + DirectAdd(Separators) - Selected - Cuts(Constraints)
 *       Selected = Applied(Separators) + Applied(Constraints)
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintCutselectorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintCutselectorStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Cutselectors       :   ExecTime  SetupTime      Calls  RootCalls   Selected     Forced   Filtered  RootSelec   RootForc   RootFilt \n");

   /* sort cutsels w.r.t. their priority */
   SCIPsetSortCutsels(scip->set);

   for( i = 0; i < scip->set->ncutsels; ++i )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
         SCIPcutselGetName(scip->set->cutsels[i]),
         SCIPcutselGetTime(scip->set->cutsels[i]),
         SCIPcutselGetSetupTime(scip->set->cutsels[i]),
         SCIPcutselGetNCalls(scip->set->cutsels[i]),
         SCIPcutselGetNRootCalls(scip->set->cutsels[i]),
         SCIPcutselGetNRootCuts(scip->set->cutsels[i]) + SCIPcutselGetNLocalCuts(scip->set->cutsels[i]),
         SCIPcutselGetNRootForcedCuts(scip->set->cutsels[i]) + SCIPcutselGetNLocalForcedCuts(scip->set->cutsels[i]),
         SCIPcutselGetNRootCutsFiltered(scip->set->cutsels[i]) + SCIPcutselGetNLocalCutsFiltered(scip->set->cutsels[i]),
         SCIPcutselGetNRootCuts(scip->set->cutsels[i]),
         SCIPcutselGetNRootForcedCuts(scip->set->cutsels[i]),
         SCIPcutselGetNRootCutsFiltered(scip->set->cutsels[i])
         );
   }
}

/** collects cutselector statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectCutselectorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* cutselectors;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectCutselectorStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Create a subtree for cutselectors */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &cutselectors, "plugins", scip->set->ncutsels) );

   /* Sort cutselectors by priority */
   SCIPsetSortCutsels(scip->set);

   /* Collect statistics for each cutselector */
   for( i = 0; i < scip->set->ncutsels; ++i )
   {
      SCIP_CUTSEL* cutsel = scip->set->cutsels[i];
      SCIP_DATATREE* cutseldata;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, cutselectors, &cutseldata, SCIPcutselGetName(cutsel), 11) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, cutseldata, "description", SCIPcutselGetDesc(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, cutseldata, "exec_time", SCIPcutselGetTime(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, cutseldata, "setup_time", SCIPcutselGetSetupTime(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "calls", SCIPcutselGetNCalls(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "root_calls", SCIPcutselGetNRootCalls(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "selected", SCIPcutselGetNRootCuts(cutsel) + SCIPcutselGetNLocalCuts(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "forced", SCIPcutselGetNRootForcedCuts(cutsel) + SCIPcutselGetNLocalForcedCuts(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "filtered", SCIPcutselGetNRootCutsFiltered(cutsel) + SCIPcutselGetNLocalCutsFiltered(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "root_selected", SCIPcutselGetNRootCuts(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "root_forced", SCIPcutselGetNRootForcedCuts(cutsel)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, cutseldata, "root_filtered", SCIPcutselGetNRootCutsFiltered(cutsel)) );
   }

   return SCIP_OKAY;
}

/** outputs pricer statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintPricerStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintPricerStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Pricers            :   ExecTime  SetupTime      Calls       Vars\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  problem variables: %10.2f          - %10d %10d\n",
      SCIPpricestoreGetProbPricingTime(scip->pricestore),
      SCIPpricestoreGetNProbPricings(scip->pricestore),
      SCIPpricestoreGetNProbvarsFound(scip->pricestore));

   /* sort pricers w.r.t. their name */
   SCIPsetSortPricersName(scip->set);

   for( i = 0; i < scip->set->nactivepricers; ++i )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10.2f %10d %10d\n",
         SCIPpricerGetName(scip->set->pricers[i]),
         SCIPpricerGetTime(scip->set->pricers[i]),
         SCIPpricerGetSetupTime(scip->set->pricers[i]),
         SCIPpricerGetNCalls(scip->set->pricers[i]),
         SCIPpricerGetNVarsFound(scip->set->pricers[i]));
   }
}

/** collects pricer statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectPricerStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* pricers;
   SCIP_DATATREE* probvars;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectPricerStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Create a subtree for pricers */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &pricers, "pricers", 1 + scip->set->nactivepricers) );

   /* Collect statistics for problem variables */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, pricers, &probvars, "problem_variables", 3) );

   SCIP_CALL( SCIPinsertDatatreeReal(scip, probvars, "exec_time", SCIPpricestoreGetProbPricingTime(scip->pricestore)) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, probvars, "calls", SCIPpricestoreGetNProbPricings(scip->pricestore)) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, probvars, "vars_found", SCIPpricestoreGetNProbvarsFound(scip->pricestore)) );

   /* Sort pricers by name */
   SCIPsetSortPricersName(scip->set);

   /* Collect statistics for each active pricer */
   for( i = 0; i < scip->set->nactivepricers; ++i )
   {
      SCIP_PRICER* pricer = scip->set->pricers[i];
      SCIP_DATATREE* pricerdata;

      SCIP_CALL( SCIPcreateDatatreeInTree( scip, pricers, &pricerdata, SCIPpricerGetName(pricer), 4) );

      SCIP_CALL( SCIPinsertDatatreeReal(scip, pricerdata, "exec_time", SCIPpricerGetTime(pricer)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, pricerdata, "setup_time", SCIPpricerGetSetupTime(pricer)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, pricerdata, "calls", SCIPpricerGetNCalls(pricer)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, pricerdata, "vars_found", SCIPpricerGetNVarsFound(pricer)) );
   }

   return SCIP_OKAY;
}

/** outputs branching rule statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintBranchruleStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintBranchruleStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Branching Rules    :   ExecTime  SetupTime   BranchLP  BranchExt   BranchPS    Cutoffs    DomReds       Cuts      Conss   Children\n");

   /* sort branching rules  w.r.t. their name */
   SCIPsetSortBranchrulesName(scip->set);

   for( i = 0; i < scip->set->nbranchrules; ++i )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
         SCIPbranchruleGetName(scip->set->branchrules[i]),
         SCIPbranchruleGetTime(scip->set->branchrules[i]),
         SCIPbranchruleGetSetupTime(scip->set->branchrules[i]),
         SCIPbranchruleGetNLPCalls(scip->set->branchrules[i]),
         SCIPbranchruleGetNExternCalls(scip->set->branchrules[i]),
         SCIPbranchruleGetNPseudoCalls(scip->set->branchrules[i]),
         SCIPbranchruleGetNCutoffs(scip->set->branchrules[i]),
         SCIPbranchruleGetNDomredsFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNCutsFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNConssFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNChildren(scip->set->branchrules[i]));
   }
}

/** collects branching rule statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectBranchruleStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* branchruletree;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(datatree != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPcollectBranchruleStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &branchruletree, "plugins", scip->set->nbranchrules) );

   /* sort branching rules  w.r.t. their name */
   SCIPsetSortBranchrulesName(scip->set);

   for( i = 0; i < scip->set->nbranchrules; ++i )
   {
      SCIP_DATATREE* branchrule = NULL;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, branchruletree, &branchrule, SCIPbranchruleGetName(scip->set->branchrules[i]), 11) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, branchrule, "description", SCIPbranchruleGetDesc(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, branchrule, "time", SCIPbranchruleGetTime(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, branchrule, "setuptime", SCIPbranchruleGetSetupTime(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "nlpcalls", SCIPbranchruleGetNLPCalls(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "nexterncalls", SCIPbranchruleGetNExternCalls(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "npscalls", SCIPbranchruleGetNPseudoCalls(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "ncutoffs", SCIPbranchruleGetNCutoffs(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "ndomreds", SCIPbranchruleGetNDomredsFound(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "ncutsfound", SCIPbranchruleGetNCutsFound(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "nconssfound", SCIPbranchruleGetNConssFound(scip->set->branchrules[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, branchrule, "nchildren", SCIPbranchruleGetNChildren(scip->set->branchrules[i])) );
   }

   return SCIP_OKAY;
}

/** outputs heuristics statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintHeuristicStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int ndivesets = 0;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->tree != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintHeuristicStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Primal Heuristics  :   ExecTime  SetupTime      Calls      Found       Best\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  LP solutions     : %10.2f          -          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
      SCIPclockGetTime(scip->stat->lpsoltime),
      scip->stat->nlpsolsfound, scip->stat->nlpbestsolsfound);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  relax solutions  : %10.2f          -          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
      SCIPclockGetTime(scip->stat->relaxsoltime),
      scip->stat->nrelaxsolsfound, scip->stat->nrelaxbestsolsfound);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  pseudo solutions : %10.2f          -          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
      SCIPclockGetTime(scip->stat->pseudosoltime),
      scip->stat->npssolsfound, scip->stat->npsbestsolsfound);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  strong branching : %10.2f          -          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
      SCIPclockGetTime(scip->stat->sbsoltime),
      scip->stat->nsbsolsfound, scip->stat->nsbbestsolsfound);

   /* sort heuristics w.r.t. their names */
   SCIPsetSortHeursName(scip->set);

   for( i = 0; i < scip->set->nheurs; ++i )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
         SCIPheurGetName(scip->set->heurs[i]),
         SCIPheurGetTime(scip->set->heurs[i]),
         SCIPheurGetSetupTime(scip->set->heurs[i]),
         SCIPheurGetNCalls(scip->set->heurs[i]),
         SCIPheurGetNSolsFound(scip->set->heurs[i]),
         SCIPheurGetNBestSolsFound(scip->set->heurs[i]));

      /* count heuristics that use diving; needed to determine output later */
      ndivesets += SCIPheurGetNDivesets(scip->set->heurs[i]);
   }

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  other solutions  :          -          -          - %10" SCIP_LONGINT_FORMAT "          -\n",
      scip->stat->nexternalsolsfound);

   if ( ndivesets > 0 && scip->set->misc_showdivingstats )
   {
      int c;
      SCIP_DIVECONTEXT divecontexts[] = {SCIP_DIVECONTEXT_SINGLE, SCIP_DIVECONTEXT_ADAPTIVE, SCIP_DIVECONTEXT_SCHEDULER};

      /* print statistics for all three contexts individually */
      for( c = 0; c < 3; ++c )
      {
         SCIP_DIVECONTEXT divecontext = divecontexts[c];

         if( divecontext == SCIP_DIVECONTEXT_SINGLE )
         {
            SCIPmessageFPrintInfo(scip->messagehdlr, file,
               "Diving %-12s:      Calls      Nodes   LP Iters Backtracks  Conflicts   MinDepth   MaxDepth   AvgDepth  RoundSols  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt\n", "(single)");
         }
         else
         {
            SCIPmessageFPrintInfo(scip->messagehdlr, file,
               "Diving %-12s:      Calls      Nodes   LP Iters Backtracks  Conflicts   MinDepth   MaxDepth   AvgDepth  RoundSols  NLeafSols  MinSolDpt  MaxSolDpt  AvgSolDpt\n",
               divecontext == SCIP_DIVECONTEXT_ADAPTIVE ? "(adaptive)" : "(scheduler)");
         }

         for( i = 0; i < scip->set->nheurs; ++i )
         {
            int s;
            for( s = 0; s < SCIPheurGetNDivesets(scip->set->heurs[i]); ++s )
            {
               SCIP_DIVESET* diveset = SCIPheurGetDivesets(scip->set->heurs[i])[s];

               SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10d",
                        SCIPdivesetGetName(diveset),
                        SCIPdivesetGetNCalls(diveset, divecontext));
               if( SCIPdivesetGetNCalls(diveset, divecontext) > 0 )
               {
                  SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10d %10d %10.1f %10" SCIP_LONGINT_FORMAT,
                           SCIPdivesetGetNProbingNodes(diveset, divecontext),
                           SCIPdivesetGetNLPIterations(diveset, divecontext),
                           SCIPdivesetGetNBacktracks(diveset, divecontext),
                           SCIPdivesetGetNConflicts(diveset, divecontext),
                           SCIPdivesetGetMinDepth(diveset, divecontext),
                           SCIPdivesetGetMaxDepth(diveset, divecontext),
                           SCIPdivesetGetAvgDepth(diveset, divecontext),
                           SCIPdivesetGetNSols(diveset, divecontext) - SCIPdivesetGetNSolutionCalls(diveset, divecontext));

                  if( SCIPdivesetGetNSolutionCalls(diveset, divecontext) > 0 )
                  {
                     SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10d %10d %10d %10.1f\n",
                              SCIPdivesetGetNSolutionCalls(diveset, divecontext),
                              SCIPdivesetGetMinSolutionDepth(diveset, divecontext),
                              SCIPdivesetGetMaxSolutionDepth(diveset, divecontext),
                              SCIPdivesetGetAvgSolutionDepth(diveset, divecontext));
                  }
                  else
                     SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -          -          -          -\n");
               }
               else
                  SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -          -          -          -          -          -          -          -          -          -          -          -\n");
            }
         }
      }
   }
}

/** collects heuristics statistics into SCIP_DATATREE
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcollectHeuristicStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* heursstore;
   SCIP_DATATREE* lpstore;
   SCIP_DATATREE* relaxstore;
   SCIP_DATATREE* pseudostore;
   SCIP_DATATREE* sbstore;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->tree != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectHeuristicStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* create the heuristics datatree */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &heursstore, "plugins", scip->set->nheurs + 5) );

   /* add LP solutions statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, heursstore, &lpstore, "lp_solutions", 4) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, lpstore, "time", SCIPclockGetTime(scip->stat->lpsoltime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, lpstore, "solutions_found", scip->stat->nlpsolsfound) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, lpstore, "best_solutions_found", scip->stat->nlpbestsolsfound) );

   /* add relaxation solutions statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, heursstore, &relaxstore, "relax_solutions", 4) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, relaxstore, "time", SCIPclockGetTime(scip->stat->relaxsoltime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxstore, "solutions_found", scip->stat->nrelaxsolsfound) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxstore, "best_solutions_found", scip->stat->nrelaxbestsolsfound) );

   /* add pseudo solutions statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, heursstore, &pseudostore, "pseudo_solutions", 4) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, pseudostore, "time", SCIPclockGetTime(scip->stat->pseudosoltime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, pseudostore, "solutions_found", scip->stat->npssolsfound) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, pseudostore, "best_solutions_found", scip->stat->npsbestsolsfound) );

   /* add strong branching statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, heursstore, &sbstore, "strong_branching", 4) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, sbstore, "time", SCIPclockGetTime(scip->stat->sbsoltime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, sbstore, "solutions_found", scip->stat->nsbsolsfound) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, sbstore, "best_solutions_found", scip->stat->nsbbestsolsfound) );

   /* Sort heuristics by name */
   SCIPsetSortHeursName(scip->set);

   for( int i = 0; i < scip->set->nheurs; ++i )
   {
      SCIP_DATATREE* heurstore;
      SCIP_CALL( SCIPcreateDatatreeInTree(scip, heursstore, &heurstore, SCIPheurGetName(scip->set->heurs[i]), 6) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, heurstore, "description", SCIPheurGetDesc(scip->set->heurs[i])) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, heurstore, "time", SCIPheurGetTime(scip->set->heurs[i])) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, heurstore, "setup_time", SCIPheurGetSetupTime(scip->set->heurs[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, heurstore, "calls", SCIPheurGetNCalls(scip->set->heurs[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, heurstore, "solutions_found", SCIPheurGetNSolsFound(scip->set->heurs[i])) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, heurstore, "best_solutions_found", SCIPheurGetNBestSolsFound(scip->set->heurs[i])) );
   }

   SCIP_DATATREE* othersolutions;
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, heursstore, &othersolutions, "other_solutions", 1) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, othersolutions, "solutions_found",  scip->stat->nexternalsolsfound) );

   /* collect diving statistics if applicable */
   if( scip->set->misc_showdivingstats )
   {
      SCIP_DATATREE* divingstore;
      SCIP_DIVECONTEXT divecontexts[] = { SCIP_DIVECONTEXT_SINGLE, SCIP_DIVECONTEXT_ADAPTIVE, SCIP_DIVECONTEXT_SCHEDULER };

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &divingstore, "diving_statistics", 3) );

      for( int c = 0; c < 3; ++c )
      {
         SCIP_DIVECONTEXT divecontext = divecontexts[c];
         SCIP_DATATREE* contextstore;

         SCIP_CALL( SCIPcreateDatatreeInTree(scip, divingstore, &contextstore,
            divecontext == SCIP_DIVECONTEXT_SINGLE ? "single" : divecontext == SCIP_DIVECONTEXT_ADAPTIVE ? "adaptive" : "scheduler",
            scip->set->nheurs) );

         for( int i = 0; i < scip->set->nheurs; ++i )
         {
            int ndivesets = SCIPheurGetNDivesets(scip->set->heurs[i]);
            for( int s = 0; s < ndivesets; ++s )
            {
               SCIP_DATATREE* divesetstore;
               SCIP_DIVESET* diveset = SCIPheurGetDivesets(scip->set->heurs[i])[s];

               SCIP_CALL( SCIPcreateDatatreeInTree(scip, contextstore, &divesetstore, SCIPdivesetGetName(diveset), 5) );

               SCIP_CALL( SCIPinsertDatatreeInt(scip, divesetstore, "calls", SCIPdivesetGetNCalls(diveset, divecontext)) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, divesetstore, "nodes", SCIPdivesetGetNProbingNodes(diveset, divecontext)) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, divesetstore, "lp_iters", SCIPdivesetGetNLPIterations(diveset, divecontext)) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, divesetstore, "backtracks", SCIPdivesetGetNBacktracks(diveset, divecontext)) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, divesetstore, "conflicts", SCIPdivesetGetNConflicts(diveset, divecontext)) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** outputs compression statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintCompressionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintCompressionStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* only print compression statistics if tree reoptimization is enabled */
   if( !scip->set->reopt_enable )
      return;

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Tree Compressions  :   ExecTime  SetupTime      Calls      Found\n");

   /* sort compressions w.r.t. their names */
   SCIPsetSortComprsName(scip->set);

   for( i = 0; i < scip->set->ncomprs; ++i )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
         SCIPcomprGetName(scip->set->comprs[i]),
         SCIPcomprGetTime(scip->set->comprs[i]),
         SCIPcomprGetSetupTime(scip->set->comprs[i]),
         SCIPcomprGetNCalls(scip->set->comprs[i]),
         SCIPcomprGetNFound(scip->set->comprs[i]));
   }
}

/** collects compression statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectCompressionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* compressions;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectCompressionStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Check if tree reoptimization is enabled */
   if( !scip->set->reopt_enable )
      return SCIP_OKAY;

   /* Create a subtree for tree compressions */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &compressions, "tree_compressions", scip->set->ncomprs) );

   /* Sort compressions by name */
   SCIPsetSortComprsName(scip->set);

   /* Collect statistics for each compression method */
   for( i = 0; i < scip->set->ncomprs; ++i )
   {
      SCIP_COMPR* compr = scip->set->comprs[i];
      SCIP_DATATREE* comprdata;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, compressions, &comprdata, SCIPcomprGetName(compr), 5) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, comprdata, "description", SCIPcomprGetDesc(compr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, comprdata, "exec_time", SCIPcomprGetTime(compr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, comprdata, "setup_time", SCIPcomprGetSetupTime(compr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, comprdata, "n_calls", SCIPcomprGetNCalls(compr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, comprdata, "n_found", SCIPcomprGetNFound(compr)) );
   }

   return SCIP_OKAY;
}

/** outputs LP statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->lp != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintLPStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It    ItLimit\n");

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  primal LP        : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->primallptime),
      scip->stat->nprimallps + scip->stat->nprimalzeroitlps,
      scip->stat->nprimallpiterations,
      scip->stat->nprimallps > 0 ? (SCIP_Real)scip->stat->nprimallpiterations/(SCIP_Real)scip->stat->nprimallps : 0.0);
   if( SCIPclockGetTime(scip->stat->primallptime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", (SCIP_Real)scip->stat->nprimallpiterations/SCIPclockGetTime(scip->stat->primallptime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10" SCIP_LONGINT_FORMAT "\n",
      scip->stat->primalzeroittime,
      scip->stat->nprimalzeroitlps);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  dual LP          : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->duallptime),
      scip->stat->nduallps + scip->stat->ndualzeroitlps,
      scip->stat->nduallpiterations,
      scip->stat->nduallps > 0 ? (SCIP_Real)scip->stat->nduallpiterations/(SCIP_Real)scip->stat->nduallps : 0.0);
   if( SCIPclockGetTime(scip->stat->duallptime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", (SCIP_Real)scip->stat->nduallpiterations/SCIPclockGetTime(scip->stat->duallptime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10" SCIP_LONGINT_FORMAT "\n",
      scip->stat->dualzeroittime,
      scip->stat->ndualzeroitlps);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  lex dual LP      : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->lexduallptime),
      scip->stat->nlexduallps,
      scip->stat->nlexduallpiterations,
      scip->stat->nlexduallps > 0 ? (SCIP_Real)scip->stat->nlexduallpiterations/(SCIP_Real)scip->stat->nlexduallps : 0.0);
   if( SCIPclockGetTime(scip->stat->lexduallptime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f\n", (SCIP_Real)scip->stat->nlexduallpiterations/SCIPclockGetTime(scip->stat->lexduallptime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -\n");

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  barrier LP       : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->barrierlptime),
      scip->stat->nbarrierlps,
      scip->stat->nbarrierlpiterations,
      scip->stat->nbarrierlps > 0 ? (SCIP_Real)scip->stat->nbarrierlpiterations/(SCIP_Real)scip->stat->nbarrierlps : 0.0);
   if( SCIPclockGetTime(scip->stat->barrierlptime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", (SCIP_Real)scip->stat->nbarrierlpiterations/SCIPclockGetTime(scip->stat->barrierlptime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f %10" SCIP_LONGINT_FORMAT "\n",
      scip->stat->barrierzeroittime,
      scip->stat->nbarrierzeroitlps);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  resolve instable : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->resolveinstablelptime),
      scip->stat->nresolveinstablelps,
      scip->stat->nresolveinstablelpiters,
      scip->stat->nresolveinstablelps > 0 ? (SCIP_Real)scip->stat->nresolveinstablelpiters/(SCIP_Real)scip->stat->nresolveinstablelps : 0.0);
   if( SCIPclockGetTime(scip->stat->resolveinstablelptime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f\n", (SCIP_Real)scip->stat->nresolveinstablelpiters/SCIPclockGetTime(scip->stat->resolveinstablelptime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -\n");

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  diving/probing LP: %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->divinglptime),
      scip->stat->ndivinglps,
      scip->stat->ndivinglpiterations,
      scip->stat->ndivinglps > 0 ? (SCIP_Real)scip->stat->ndivinglpiterations/(SCIP_Real)scip->stat->ndivinglps : 0.0);
   if( SCIPclockGetTime(scip->stat->divinglptime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f\n", (SCIP_Real)scip->stat->ndivinglpiterations/SCIPclockGetTime(scip->stat->divinglptime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -\n");

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  strong branching : %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->strongbranchtime),
      scip->stat->nstrongbranchs,
      scip->stat->nsblpiterations,
      scip->stat->nstrongbranchs > 0 ? (SCIP_Real)scip->stat->nsblpiterations/(SCIP_Real)scip->stat->nstrongbranchs : 0.0);
   if( SCIPclockGetTime(scip->stat->strongbranchtime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", (SCIP_Real)scip->stat->nsblpiterations/SCIPclockGetTime(scip->stat->strongbranchtime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -          - %10" SCIP_LONGINT_FORMAT "\n", scip->stat->nsbtimesiterlimhit);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "    (at root node) :          - %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f          -\n",
      scip->stat->nrootstrongbranchs,
      scip->stat->nrootsblpiterations,
      scip->stat->nrootstrongbranchs > 0
      ? (SCIP_Real)scip->stat->nrootsblpiterations/(SCIP_Real)scip->stat->nrootstrongbranchs : 0.0);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  conflict analysis: %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f",
      SCIPclockGetTime(scip->stat->conflictlptime),
      scip->stat->nconflictlps,
      scip->stat->nconflictlpiterations,
      scip->stat->nconflictlps > 0 ? (SCIP_Real)scip->stat->nconflictlpiterations/(SCIP_Real)scip->stat->nconflictlps : 0.0);
   if( SCIPclockGetTime(scip->stat->conflictlptime) >= 0.01 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f\n", (SCIP_Real)scip->stat->nconflictlpiterations/SCIPclockGetTime(scip->stat->conflictlptime));
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "          -\n");
}

/** collects LP statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* primal;
   SCIP_DATATREE* dual;
   SCIP_DATATREE* barrier;
   SCIP_DATATREE* lexdual;
   SCIP_DATATREE* resolve;
   SCIP_DATATREE* strongbranch;
   SCIP_DATATREE* strongbranchroot;
   SCIP_DATATREE* conflict;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->lp != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectLPStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Primal LP statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &primal, "primal_lp", 4) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, primal, "time", SCIPclockGetTime(scip->stat->primallptime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, primal, "calls", scip->stat->nprimallps + scip->stat->nprimalzeroitlps) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, primal, "iterations", scip->stat->nprimallpiterations) );
   if( scip->stat->nprimallps > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, primal, "iter_per_call", (SCIP_Real)scip->stat->nprimallpiterations / scip->stat->nprimallps) );
   }

   /* Dual LP statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &dual, "dual_lp", 4) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, dual, "time", SCIPclockGetTime(scip->stat->duallptime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, dual, "calls", scip->stat->nduallps + scip->stat->ndualzeroitlps) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, dual, "iterations", scip->stat->nduallpiterations) );
   if( scip->stat->nduallps > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, dual, "iter_per_call", (SCIP_Real)scip->stat->nduallpiterations / scip->stat->nduallps) );
   }

   /* Barrier LP statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &barrier, "barrier_lp", 4) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, barrier, "time", SCIPclockGetTime(scip->stat->barrierlptime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, barrier, "calls", scip->stat->nbarrierlps) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, barrier, "iterations", scip->stat->nbarrierlpiterations) );
   if( scip->stat->nbarrierlps > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, barrier, "iter_per_call", (SCIP_Real)scip->stat->nbarrierlpiterations / scip->stat->nbarrierlps) );
   }

   /* Lex dual LP statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &lexdual, "lex_dual_lp", 5) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, lexdual, "time", SCIPclockGetTime(scip->stat->lexduallptime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, lexdual, "calls", scip->stat->nlexduallps) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, lexdual, "iterations", scip->stat->nlexduallpiterations) );
   if( scip->stat->nlexduallps > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, lexdual, "iter_per_call", (SCIP_Real)scip->stat->nlexduallpiterations / scip->stat->nlexduallps) );
   }
   if( SCIPclockGetTime(scip->stat->lexduallptime) >= 0.01 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, lexdual, "iter_per_time", (SCIP_Real)scip->stat->nlexduallpiterations / SCIPclockGetTime(scip->stat->lexduallptime)) );
   }

   /* Resolving LP statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &resolve, "resolve_lp", 5) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, resolve, "time", SCIPclockGetTime(scip->stat->resolveinstablelptime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, resolve, "calls", scip->stat->nresolveinstablelps) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, resolve, "iterations", scip->stat->nresolveinstablelpiters) );
   if( scip->stat->nresolveinstablelps > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, resolve, "iter_per_call", (SCIP_Real)scip->stat->nresolveinstablelpiters / scip->stat->nresolveinstablelps) );
   }
   if( SCIPclockGetTime(scip->stat->resolveinstablelptime) >= 0.01 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, resolve, "iter_per_time", (SCIP_Real)scip->stat->nresolveinstablelpiters / SCIPclockGetTime(scip->stat->resolveinstablelptime)) );
   }

   /* Strong branching LP statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &strongbranch, "strongbranch_lp", 5) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, strongbranch, "time", SCIPclockGetTime(scip->stat->strongbranchtime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, strongbranch, "calls", scip->stat->nstrongbranchs) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, strongbranch, "iterations", scip->stat->nsblpiterations) );
   if( scip->stat->nstrongbranchs > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, strongbranch, "iter_per_call", (SCIP_Real)scip->stat->nsblpiterations / scip->stat->nstrongbranchs) );
   }
   if( SCIPclockGetTime(scip->stat->strongbranchtime) >= 0.01 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, strongbranch, "iter_per_time", (SCIP_Real)scip->stat->nsblpiterations / SCIPclockGetTime(scip->stat->strongbranchtime)) );
   }

   /* Strong branching at root node statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &strongbranchroot, "strongbranch_root_lp", 3) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, strongbranchroot, "calls", scip->stat->nrootstrongbranchs) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, strongbranchroot, "iterations", scip->stat->nrootsblpiterations) );
   if( scip->stat->nrootstrongbranchs > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, strongbranchroot, "iter_per_call", (SCIP_Real)scip->stat->nrootsblpiterations / scip->stat->nrootstrongbranchs) );
   }

   /* Conflict analysis LP statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &conflict, "conflict_lp", 5) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, conflict, "time", SCIPclockGetTime(scip->stat->conflictlptime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, conflict, "calls", scip->stat->nconflictlps) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, conflict, "iterations", scip->stat->nconflictlpiterations) );
   if( scip->stat->nconflictlps > 0 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conflict, "iter_per_call", (SCIP_Real)scip->stat->nconflictlpiterations / scip->stat->nconflictlps) );
   }
   if( SCIPclockGetTime(scip->stat->conflictlptime) >= 0.01 )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, conflict, "iter_per_time", (SCIP_Real)scip->stat->nconflictlpiterations / SCIPclockGetTime(scip->stat->conflictlptime)) );
   }

   return SCIP_OKAY;
}

/** outputs NLP statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintNLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int nnlrowlinear;
   int nnlrowconvexineq;
   int nnlrownonconvexineq;
   int nnlrownonlineareq;

   assert(scip != NULL);
   assert(scip->stat != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintNLPStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
      return;

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "NLP relaxation     :\n");

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  solve time       : %10.2f (%" SCIP_LONGINT_FORMAT " calls)\n",
      SCIPclockGetTime(scip->stat->nlpsoltime),
      scip->stat->nnlps);

   SCIP_CALL_ABORT( SCIPgetNLPNlRowsStat(scip, &nnlrowlinear, &nnlrowconvexineq, &nnlrownonconvexineq, &nnlrownonlineareq) );
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  convexity        : %10s (%d linear rows, %d convex ineq., %d nonconvex ineq., %d nonlinear eq. or two-sided ineq.)\n",
      (nnlrownonconvexineq == 0 && nnlrownonlineareq == 0) ? "convex" : "nonconvex",
      nnlrowlinear, nnlrowconvexineq, nnlrownonconvexineq, nnlrownonlineareq);
}

/** collects NLP statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectNLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   int nnlrowlinear;
   int nnlrowconvexineq;
   int nnlrownonconvexineq;
   int nnlrownonlineareq;

   assert(scip != NULL);
   assert(datatree != NULL);
   assert(scip->stat != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectNLPStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "solve_time", SCIPclockGetTime(scip->stat->nlpsoltime)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "calls", scip->stat->nnlps) );

   /* Get convexity and row statistics */
   SCIP_CALL( SCIPgetNLPNlRowsStat(scip, &nnlrowlinear, &nnlrowconvexineq, &nnlrownonconvexineq, &nnlrownonlineareq) );

   SCIP_CALL( SCIPinsertDatatreeString(scip, datatree, "convexity",
      (nnlrownonconvexineq == 0 && nnlrownonlineareq == 0) ? "convex" : "nonconvex") );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "linear_rows", nnlrowlinear) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "convex_ineq", nnlrowconvexineq) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "nonconvex_ineq", nnlrownonconvexineq) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "nonlinear_eq", nnlrownonlineareq) );

   return SCIP_OKAY;
}

/** outputs relaxator statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintRelaxatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintRelaxatorStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->set->nrelaxs == 0 )
      return;

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Relaxators         :       Time      Calls    Cutoffs ImprBounds   ImprTime ReducedDom  Separated AddedConss\n");

   /* sort relaxators w.r.t. their name */
   SCIPsetSortRelaxsName(scip->set);

   for( i = 0; i < scip->set->nrelaxs; ++i )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
         SCIPrelaxGetName(scip->set->relaxs[i]),
         SCIPrelaxGetTime(scip->set->relaxs[i]),
         SCIPrelaxGetNCalls(scip->set->relaxs[i]),
         SCIPrelaxGetNCutoffs(scip->set->relaxs[i]),
         SCIPrelaxGetNImprovedLowerbound(scip->set->relaxs[i]),
         SCIPrelaxGetImprovedLowerboundTime(scip->set->relaxs[i]),
         SCIPrelaxGetNReducedDomains(scip->set->relaxs[i]),
         SCIPrelaxGetNSeparatedCuts(scip->set->relaxs[i]),
         SCIPrelaxGetNAddedConss(scip->set->relaxs[i])
         );
   }
}

/** collects relaxator statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectRelaxatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* relaxators;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectRelaxatorStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->set->nrelaxs == 0 )
      return SCIP_OKAY;

   /* Create a subtree for relaxators */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &relaxators, "relaxators", scip->set->nrelaxs) );

   /* Sort relaxators by name */
   SCIPsetSortRelaxsName(scip->set);

   /* Collect statistics for each relaxator */
   for( i = 0; i < scip->set->nrelaxs; ++i )
   {
      SCIP_RELAX* relax = scip->set->relaxs[i];
      SCIP_DATATREE* relaxdata;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, relaxators, &relaxdata, SCIPrelaxGetName(relax), 8) );

      SCIP_CALL( SCIPinsertDatatreeReal(scip, relaxdata, "time", SCIPrelaxGetTime(relax)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxdata, "calls", SCIPrelaxGetNCalls(relax)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxdata, "cutoffs", SCIPrelaxGetNCutoffs(relax)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxdata, "improved_lowerbound", SCIPrelaxGetNImprovedLowerbound(relax)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, relaxdata, "improved_lowerbound_time", SCIPrelaxGetImprovedLowerboundTime(relax)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxdata, "reduced_domains", SCIPrelaxGetNReducedDomains(relax)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxdata, "separated_cuts", SCIPrelaxGetNSeparatedCuts(relax)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, relaxdata, "added_constraints", SCIPrelaxGetNAddedConss(relax)) );
   }

   return SCIP_OKAY;
}

/** outputs tree statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintTreeStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->tree != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintTreeStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "B&B Tree           :\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  number of runs   : %10d\n", scip->stat->nruns);
   SCIPmessageFPrintInfo(scip->messagehdlr, file,
      "  nodes            : %10" SCIP_LONGINT_FORMAT " (%" SCIP_LONGINT_FORMAT " internal, %" SCIP_LONGINT_FORMAT " leaves)\n",
      scip->stat->nnodes, scip->stat->ninternalnodes, scip->stat->nnodes - scip->stat->ninternalnodes );
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  feasible leaves  : %10" SCIP_LONGINT_FORMAT "\n", scip->stat->nfeasleaves);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  infeas. leaves   : %10" SCIP_LONGINT_FORMAT "\n", scip->stat->ninfeasleaves);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  objective leaves : %10" SCIP_LONGINT_FORMAT "\n", scip->stat->nobjleaves);
   SCIPmessageFPrintInfo(scip->messagehdlr, file,
      "  nodes (total)    : %10" SCIP_LONGINT_FORMAT " (%" SCIP_LONGINT_FORMAT " internal, %" SCIP_LONGINT_FORMAT " leaves)\n",
      scip->stat->ntotalnodes, scip->stat->ntotalinternalnodes, scip->stat->ntotalnodes - scip->stat->ntotalinternalnodes);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  nodes left       : %10d\n", SCIPtreeGetNNodes(scip->tree));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  max depth        : %10d\n", scip->stat->maxdepth);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  max depth (total): %10d\n", scip->stat->maxtotaldepth);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  backtracks       : %10" SCIP_LONGINT_FORMAT " (%.1f%%)\n", scip->stat->nbacktracks,
      scip->stat->nnodes > 0 ? 100.0 * (SCIP_Real)scip->stat->nbacktracks / (SCIP_Real)scip->stat->nnodes : 0.0);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  early backtracks : %10" SCIP_LONGINT_FORMAT " (%.1f%%)\n", scip->stat->nearlybacktracks,
       scip->stat->nbacktracks > 0 ? 100.0 * (SCIP_Real)scip->stat->nearlybacktracks / (SCIP_Real)scip->stat->nbacktracks : 0.0);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  nodes exc. ref.  : %10" SCIP_LONGINT_FORMAT " (%.1f%%)\n", scip->stat->nnodesaboverefbound,
       scip->stat->nnodes > 0 ? 100.0 * (SCIP_Real)scip->stat->nnodesaboverefbound / (SCIP_Real)scip->stat->nnodes : 0.0);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  delayed cutoffs  : %10" SCIP_LONGINT_FORMAT "\n", scip->stat->ndelayedcutoffs);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  repropagations   : %10" SCIP_LONGINT_FORMAT " (%" SCIP_LONGINT_FORMAT " domain reductions, %" SCIP_LONGINT_FORMAT " cutoffs)\n",
      scip->stat->nreprops, scip->stat->nrepropboundchgs, scip->stat->nrepropcutoffs);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  avg switch length: %10.2f\n",
      scip->stat->nnodes > 0
      ? (SCIP_Real)(scip->stat->nactivatednodes + scip->stat->ndeactivatednodes) / (SCIP_Real)scip->stat->nnodes : 0.0);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  switching time   : %10.2f\n", SCIPclockGetTime(scip->stat->nodeactivationtime));
}

/** collects tree statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectTreeStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* nodes;
   SCIP_DATATREE* totalnodes;
   SCIP_DATATREE* depth;
   SCIP_DATATREE* backtracks;
   SCIP_DATATREE* reprop;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->tree != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectTreeStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* General node statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &nodes, "nodes", 6) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, nodes, "total", scip->stat->nnodes) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, nodes, "internal", scip->stat->ninternalnodes) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, nodes, "leaves", scip->stat->nnodes - scip->stat->ninternalnodes) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, nodes, "feasible_leaves", scip->stat->nfeasleaves) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, nodes, "infeasible_leaves", scip->stat->ninfeasleaves) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, nodes, "objective_leaves", scip->stat->nobjleaves) );

   /* Total node statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &totalnodes, "total_nodes", 3) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, totalnodes, "total", scip->stat->ntotalnodes) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, totalnodes, "internal", scip->stat->ntotalinternalnodes) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, totalnodes, "leaves", scip->stat->ntotalnodes - scip->stat->ntotalinternalnodes) );

   /* Depth statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &depth, "depth", 2) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, depth, "max_depth", scip->stat->maxdepth) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, depth, "max_depth_total", scip->stat->maxtotaldepth) );

   /* Backtracks statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &backtracks, "backtracks", 4) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, backtracks, "total", scip->stat->nbacktracks) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, backtracks, "percent",
      scip->stat->nnodes > 0 ? 100.0 * (SCIP_Real)scip->stat->nbacktracks / (SCIP_Real)scip->stat->nnodes : 0.0) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, backtracks, "early", scip->stat->nearlybacktracks) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, backtracks, "early_percent",
      scip->stat->nbacktracks > 0 ? 100.0 * (SCIP_Real)scip->stat->nearlybacktracks / (SCIP_Real)scip->stat->nbacktracks : 0.0) );

   /* Repropagations statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &reprop, "repropagations", 3) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, reprop, "total", scip->stat->nreprops) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, reprop, "domain_reductions", scip->stat->nrepropboundchgs) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, reprop, "cutoffs", scip->stat->nrepropcutoffs) );

   /* Additional statistics */
   SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "nodes_left", SCIPtreeGetNNodes(scip->tree)) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "delayed_cutoffs", scip->stat->ndelayedcutoffs) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "average_switch_length",
      scip->stat->nnodes > 0 ? (SCIP_Real)(scip->stat->nactivatednodes + scip->stat->ndeactivatednodes) / (SCIP_Real)scip->stat->nnodes : 0.0) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "switching_time", SCIPclockGetTime(scip->stat->nodeactivationtime)) );

   return SCIP_OKAY;
}

/** outputs solution statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintSolutionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_Real primalbound;
   SCIP_Real dualbound;
   SCIP_Real gap;
   SCIP_Real firstprimalbound;
   SCIP_Bool objlimitreached;
   char limsolstring[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->primal != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintSolutionStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);
   gap = SCIPgetGap(scip);

   /* We output that the objective limit has been reached if the problem has been solved, no solution respecting the
    * objective limit has been found (nlimsolsfound == 0) and the primal bound is finite. Note that it still might be
    * that the original problem is infeasible, even without the objective limit, i.e., we cannot be sure that we
    * actually reached the objective limit. */
   objlimitreached = FALSE;
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED && scip->primal->nlimsolsfound == 0
      && !SCIPisInfinity(scip, primalbound) && SCIPgetStatus(scip) != SCIP_STATUS_INFORUNBD )
      objlimitreached = TRUE;

   if( scip->primal->nsolsfound != scip->primal->nlimsolsfound )
      (void) SCIPsnprintf(limsolstring, SCIP_MAXSTRLEN, ", %" SCIP_LONGINT_FORMAT " respecting the objective limit", scip->primal->nlimsolsfound);
   else
      limsolstring[0] = '\0';

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Solution           :\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Solutions found  : %10" SCIP_LONGINT_FORMAT " (%" SCIP_LONGINT_FORMAT " improvements%s)\n",
      scip->primal->nsolsfound, scip->primal->nbestsolsfound, limsolstring);

   if( SCIPsetIsInfinity(scip->set, REALABS(primalbound)) )
   {
      if( scip->set->stage == SCIP_STAGE_SOLVED )
      {
         if( scip->primal->nlimsolsfound == 0 )
         {
            if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD )
            {
               assert(!objlimitreached);
               SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Primal Bound     : infeasible or unbounded\n");
            }
            else
            {
               assert(SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE);
               if( objlimitreached )
                  SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Primal Bound     : infeasible (objective limit reached)\n");
               else
                  SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Primal Bound     : infeasible\n");
            }
         }
         else
         {
            assert(!objlimitreached);
            assert(SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED);
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Primal Bound     :  unbounded\n");
         }
      }
      else
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Primal Bound     :          -\n");
   }
   else
   {
      if( scip->primal->nlimsolsfound == 0 )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Primal Bound     : %+21.14e   (objective limit)\n", primalbound);

         /* display (best) primal bound */
         if( scip->primal->nsolsfound > 0 )
         {
            SCIP_Real bestsol;
            bestsol = SCIPsolGetObj(scip->primal->sols[0], scip->set, scip->transprob, scip->origprob);
            bestsol = SCIPretransformObj(scip, bestsol);

            SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Best Solution    : %+21.14e\n", bestsol);
         }
      }
      else
      {
         /* display first primal bound line */
         firstprimalbound = scip->stat->firstprimalbound;
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  First Solution   : %+21.14e", firstprimalbound);

         SCIPmessageFPrintInfo(scip->messagehdlr, file, "   (in run %d, after %" SCIP_LONGINT_FORMAT " nodes, %.2f seconds, depth %d, found by <%s>)\n",
            scip->stat->nrunsbeforefirst,
            scip->stat->nnodesbeforefirst,
            scip->stat->firstprimaltime,
            scip->stat->firstprimaldepth,
            ( scip->stat->firstprimalheur != NULL )
            ? ( SCIPheurGetName(scip->stat->firstprimalheur) )
            : (( scip->stat->nrunsbeforefirst == 0 ) ? "initial" : "relaxation"));

         if( SCIPisInfinity(scip, scip->stat->firstsolgap) )
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Gap First Sol.   :   infinite\n");
         else
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Gap First Sol.   : %10.2f %%\n", 100.0 * scip->stat->firstsolgap);

         if( SCIPisInfinity(scip, scip->stat->lastsolgap) )
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Gap Last Sol.    :   infinite\n");
         else
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Gap Last Sol.    : %10.2f %%\n",  100.0 * scip->stat->lastsolgap);

         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Primal Bound     : %+21.14e", primalbound);

         SCIPmessageFPrintInfo(scip->messagehdlr, file, "   (in run %d, after %" SCIP_LONGINT_FORMAT " nodes, %.2f seconds, depth %d, found by <%s>)\n",
            SCIPsolGetRunnum(scip->primal->sols[0]),
            SCIPsolGetNodenum(scip->primal->sols[0]),
            SCIPsolGetTime(scip->primal->sols[0]),
            SCIPsolGetDepth(scip->primal->sols[0]),
            SCIPsolGetHeur(scip->primal->sols[0]) != NULL
            ? SCIPheurGetName(SCIPsolGetHeur(scip->primal->sols[0]))
            : (SCIPsolGetRunnum(scip->primal->sols[0]) == 0 ? "initial" : "relaxation"));
      }
   }

   if( SCIPsetIsInfinity(scip->set, REALABS(dualbound)) )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Dual Bound       :          -\n");
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Dual Bound       : %+21.14e\n", dualbound);

   if( SCIPsetIsInfinity(scip->set, gap) )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Gap              :   infinite\n");
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Gap              : %10.2f %%\n", 100.0 * gap);

   if( scip->set->misc_calcintegral )
   {
      int s;
      const char* names[] = {
         "primal-dual",
         "primal-ref",
         "dual-ref"
      };
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "Integrals          :      Total       Avg%%\n");
      if( SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE && ! objlimitreached )
      {
        for( s = 0; s < 3; ++s )
        {
           SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17s: %10s %10s (problem infeasible)\n",
              names[s], "-", "-");
        }
      }
      else
      {
         SCIP_Real integrals[3];
         SCIP_Real solvingtime = SCIPgetSolvingTime(scip);

         if( !SCIPisFeasZero(scip, solvingtime) )
         {
            integrals[0] =  SCIPstatGetPrimalDualIntegral(scip->stat, scip->set, scip->transprob, scip->origprob, TRUE);

            if( scip->set->misc_referencevalue != SCIP_INVALID ) /*lint !e777*/
            {
               integrals[1] = SCIPstatGetPrimalReferenceIntegral(scip->stat, scip->set, scip->transprob, scip->origprob, FALSE);
               integrals[2] = SCIPstatGetDualReferenceIntegral(scip->stat, scip->set, scip->transprob, scip->origprob, FALSE);
            }
            else
               integrals[1] = integrals[2] = SCIP_INVALID;
         }
         else
         {
            BMSclearMemoryArray(integrals, 3);
         }

         /* print integrals, if computed */
         for( s = 0; s < 3; ++s )
         {
            if( integrals[s] == SCIP_INVALID ) /*lint !e777*/
               SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17s:          -          - (not evaluated)\n", names[s]);
            else
            {
               SCIP_Real avg = integrals[s] / MAX(solvingtime,1e-6);

               /* caution: this assert is non-deterministic since it depends on the solving time */
               assert(0.0 <= avg && SCIPisLE(scip, avg, 100.0));
               SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17s: %10.2f %10.2f\n", names[s], integrals[s], avg);
            }
         }
      }
   }
}

/** collects solution statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectSolutionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_Real firstprimalbound;
   SCIP_Real primalbound;
   SCIP_Real dualbound;
   SCIP_Real gap;
   SCIP_Bool objlimitreached;

   assert(scip != NULL);
   assert(datatree != NULL);
   assert(scip->stat != NULL);
   assert(scip->primal != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectSolutionStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Basic statistics */
   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);
   gap = SCIPgetGap(scip);
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "solutions_found", scip->primal->nsolsfound) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "improvements", scip->primal->nbestsolsfound) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "primal_bound", primalbound) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "dual_bound", dualbound) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "gap", gap) );

   /* Objective limit reached */
   objlimitreached = (SCIPgetStage(scip) == SCIP_STAGE_SOLVED) && (scip->primal->nlimsolsfound == 0) &&
      !SCIPisInfinity(scip, primalbound) && (SCIPgetStatus(scip) != SCIP_STATUS_INFORUNBD);
   SCIP_CALL( SCIPinsertDatatreeBool(scip, datatree, "objective_limit_reached", objlimitreached) );

   /* First solution statistics */
   firstprimalbound = scip->stat->firstprimalbound;
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "first_primal_bound", firstprimalbound) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "first_solution_time", scip->stat->firstprimaltime) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "first_solution_nodes", scip->stat->nnodesbeforefirst) );
   SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "first_solution_depth", scip->stat->firstprimaldepth) );

   return SCIP_OKAY;
}

/** outputs concurrent solver statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintConcsolverStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_CONCSOLVER** concsolvers;
   int               nconcsolvers;
   int               i;
   int               winner;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintConcsolverStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPsyncstoreIsInitialized(scip->syncstore) )
      return;

   nconcsolvers = SCIPgetNConcurrentSolvers(scip);
   concsolvers = SCIPgetConcurrentSolvers(scip);
   winner = SCIPsyncstoreGetWinner(scip->syncstore);

   if( nconcsolvers > 0 )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "Concurrent Solvers : SolvingTime    SyncTime       Nodes    LP Iters SolsShared   SolsRecvd TighterBnds TighterIntBnds\n");
      for( i = 0; i < nconcsolvers; ++i )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %c%-16s: %11.2f %11.2f %11" SCIP_LONGINT_FORMAT " %11" SCIP_LONGINT_FORMAT "%11" SCIP_LONGINT_FORMAT " %11" SCIP_LONGINT_FORMAT " %11" SCIP_LONGINT_FORMAT " %14" SCIP_LONGINT_FORMAT "\n",
            winner == i ? '*' : ' ',
            SCIPconcsolverGetName(concsolvers[i]),
            SCIPconcsolverGetSolvingTime(concsolvers[i]),
            SCIPconcsolverGetSyncTime(concsolvers[i]),
            SCIPconcsolverGetNNodes(concsolvers[i]),
            SCIPconcsolverGetNLPIterations(concsolvers[i]),
            SCIPconcsolverGetNSolsShared(concsolvers[i]),
            SCIPconcsolverGetNSolsRecvd(concsolvers[i]),
            SCIPconcsolverGetNTighterBnds(concsolvers[i]),
            SCIPconcsolverGetNTighterIntBnds(concsolvers[i])
            );
      }
   }
}

/** collects concurrent solver statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectConcsolverStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* concsolverstree;
   SCIP_CONCSOLVER** concsolvers;
   int nconcsolvers;
   int winner;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectConcsolverStatistics", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPsyncstoreIsInitialized(scip->syncstore) )
      return SCIP_OKAY;

   nconcsolvers = SCIPgetNConcurrentSolvers(scip);
   concsolvers = SCIPgetConcurrentSolvers(scip);
   winner = SCIPsyncstoreGetWinner(scip->syncstore);

   if( nconcsolvers == 0 )
      return SCIP_OKAY;

   /* Create a subtree for concurrent solvers */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &concsolverstree, "concurrent_solvers", nconcsolvers) );

   for( i = 0; i < nconcsolvers; ++i )
   {
      SCIP_CONCSOLVER* solver = concsolvers[i];
      SCIP_DATATREE* solverdata;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, concsolverstree, &solverdata, SCIPconcsolverGetName(solver), 9) );

      SCIP_CALL( SCIPinsertDatatreeInt(scip, solverdata, "is_winner", winner == i) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, solverdata, "solving_time", SCIPconcsolverGetSolvingTime(solver)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, solverdata, "sync_time", SCIPconcsolverGetSyncTime(solver)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, solverdata, "nodes", SCIPconcsolverGetNNodes(solver)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, solverdata, "lp_iterations", SCIPconcsolverGetNLPIterations(solver)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, solverdata, "solutions_shared", SCIPconcsolverGetNSolsShared(solver)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, solverdata, "solutions_received", SCIPconcsolverGetNSolsRecvd(solver)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, solverdata, "tighter_bounds", SCIPconcsolverGetNTighterBnds(solver)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, solverdata, "tighter_integer_bounds", SCIPconcsolverGetNTighterIntBnds(solver)) );
   }

   return SCIP_OKAY;
}

/** display Benders' decomposition statistics */
void SCIPprintBendersStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_BENDERS** benders;
   int nbenders;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintBendersStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPgetNActiveBenders(scip) == 0 )
      return;

   nbenders = SCIPgetNBenders(scip);
   benders = SCIPgetBenders(scip);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Benders Decomp     :   ExecTime  SetupTime      Calls      Found   Transfer   StrCalls   StrFails    StrCuts\n");
   for( i = 0; i < nbenders; ++i )
   {
      if( SCIPbendersIsActive(benders[i]) )
      {
         SCIP_BENDERSCUT** benderscuts;
         int nbenderscuts;
         int j;

         SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17.17s: %10.2f %10.2f %10d %10d %10d %10d %10d %10d\n",
            SCIPbendersGetName(scip->set->benders[i]),
            SCIPbendersGetTime(scip->set->benders[i]),
            SCIPbendersGetSetupTime(scip->set->benders[i]),
            SCIPbendersGetNCalls(scip->set->benders[i]),
            SCIPbendersGetNCutsFound(scip->set->benders[i]),
            SCIPbendersGetNTransferredCuts(scip->set->benders[i]),
            SCIPbendersGetNStrengthenCalls(scip->set->benders[i]),
            SCIPbendersGetNStrengthenFails(scip->set->benders[i]),
            SCIPbendersGetNStrengthenCutsFound(scip->set->benders[i]));

         nbenderscuts = SCIPbendersGetNBenderscuts(scip->set->benders[i]);
         benderscuts = SCIPbendersGetBenderscuts(scip->set->benders[i]);

         for( j = 0; j < nbenderscuts; j++ )
         {
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "    %-15.17s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "          -\n",
               SCIPbenderscutGetName(benderscuts[j]),
               SCIPbenderscutGetTime(benderscuts[j]),
               SCIPbenderscutGetSetupTime(benderscuts[j]),
               SCIPbenderscutGetNCalls(benderscuts[j]),
               SCIPbenderscutGetNFound(benderscuts[j]));
         }
      }
   }
}

/** collects Benders' decomposition statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectBendersStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* benderstree;
   SCIP_BENDERS** benders;
   int nbenders;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectBendersStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPgetNActiveBenders(scip) == 0 )
      return SCIP_OKAY;

   nbenders = SCIPgetNBenders(scip);
   benders = SCIPgetBenders(scip);

   /* Create a subtree for Benders statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &benderstree, "plugins", nbenders) );

   for( i = 0; i < nbenders; ++i )
   {
      if( SCIPbendersIsActive(benders[i]) )
      {
         SCIP_BENDERSCUT** benderscuts = SCIPbendersGetBenderscuts(benders[i]);
         int nbenderscuts = SCIPbendersGetNBenderscuts(benders[i]);
         SCIP_DATATREE* bendersdata;

         SCIP_CALL( SCIPcreateDatatreeInTree(scip, benderstree, &bendersdata, SCIPbendersGetName(benders[i]), 9 + nbenderscuts) );

         SCIP_CALL( SCIPinsertDatatreeString(scip, bendersdata, "description", SCIPbendersGetDesc(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, bendersdata, "exec_time", SCIPbendersGetTime(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, bendersdata, "setup_time", SCIPbendersGetSetupTime(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, bendersdata, "calls", SCIPbendersGetNCalls(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, bendersdata, "cuts_found", SCIPbendersGetNCutsFound(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, bendersdata, "cuts_transferred", SCIPbendersGetNTransferredCuts(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, bendersdata, "strength_calls", SCIPbendersGetNStrengthenCalls(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, bendersdata, "strength_failures", SCIPbendersGetNStrengthenFails(benders[i])) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, bendersdata, "strength_cuts_found", SCIPbendersGetNStrengthenCutsFound(benders[i])) );

         /* Collect statistics for Benders' cuts */
         for( int j = 0; j < nbenderscuts; ++j )
         {
            SCIP_DATATREE* benderscutdata;

            SCIP_CALL( SCIPcreateDatatreeInTree(scip, bendersdata, &benderscutdata, SCIPbenderscutGetName(benderscuts[j]), 5) );

            SCIP_CALL( SCIPinsertDatatreeString(scip, benderscutdata, "description", SCIPbenderscutGetDesc(benderscuts[j])) );
            SCIP_CALL( SCIPinsertDatatreeReal(scip, benderscutdata, "exec_time", SCIPbenderscutGetTime(benderscuts[j])) );
            SCIP_CALL( SCIPinsertDatatreeReal(scip, benderscutdata, "setup_time", SCIPbenderscutGetSetupTime(benderscuts[j])) );
            SCIP_CALL( SCIPinsertDatatreeLong(scip, benderscutdata, "calls", SCIPbenderscutGetNCalls(benderscuts[j])) );
            SCIP_CALL( SCIPinsertDatatreeLong(scip, benderscutdata, "cuts_found", SCIPbenderscutGetNFound(benderscuts[j])) );
         }
      }
   }

   return SCIP_OKAY;
}

/** outputs root statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintRootStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_Real dualboundroot;
   SCIP_Real firstdualboundroot;
   SCIP_Real firstlptime;
   SCIP_Real firstlpspeed;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->primal != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintRootStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   dualboundroot = SCIPgetDualboundRoot(scip);
   firstdualboundroot = SCIPgetFirstLPDualboundRoot(scip);
   firstlptime = SCIPgetFirstLPTime(scip);

   if( firstlptime > 0.0 )
      firstlpspeed = (SCIP_Real)scip->stat->nrootfirstlpiterations/firstlptime;
   else
      firstlpspeed = 0.0;

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Root Node          :\n");
   if( SCIPsetIsInfinity(scip->set, REALABS(firstdualboundroot)) )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  First LP value   :          -\n");
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  First LP value   : %+21.14e\n", firstdualboundroot);
   if( firstlpspeed > 0.0 )
      SCIPmessageFPrintInfo(scip->messagehdlr, file,    "  First LP Iters   : %10" SCIP_LONGINT_FORMAT " (%.2f Iter/sec)\n",
         scip->stat->nrootfirstlpiterations,
         (SCIP_Real)scip->stat->nrootfirstlpiterations/firstlptime);
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file,    "  First LP Iters   : %10" SCIP_LONGINT_FORMAT "\n", scip->stat->nrootfirstlpiterations);
   SCIPmessageFPrintInfo(scip->messagehdlr, file,    "  First LP Time    : %10.2f\n", firstlptime);

   if( SCIPsetIsInfinity(scip->set, REALABS(dualboundroot)) )
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Final Dual Bound :          -\n");
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Final Dual Bound : %+21.14e\n", dualboundroot);
   SCIPmessageFPrintInfo(scip->messagehdlr, file,    "  Final Root Iters : %10" SCIP_LONGINT_FORMAT "\n", scip->stat->nrootlpiterations);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  Root LP Estimate : ");
   if( scip->stat->rootlpbestestimate != SCIP_INVALID ) /*lint !e777*/
   {
       SCIPmessageFPrintInfo(scip->messagehdlr, file, "%+21.14e\n", SCIPretransformObj(scip, scip->stat->rootlpbestestimate));
   }
   else
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "%21s\n","-");
}

/** collects root node statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectRootStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_Real dualboundroot;
   SCIP_Real firstdualboundroot;
   SCIP_Real firstlptime;
   SCIP_Real firstlpspeed;

   assert(scip != NULL);
   assert(datatree != NULL);
   assert(scip->stat != NULL);
   assert(scip->primal != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectRootStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   dualboundroot = SCIPgetDualboundRoot(scip);
   firstdualboundroot = SCIPgetFirstLPDualboundRoot(scip);
   firstlptime = SCIPgetFirstLPTime(scip);

   if( firstlptime > 0.0 )
      firstlpspeed = (SCIP_Real)scip->stat->nrootfirstlpiterations / firstlptime;
   else
      firstlpspeed = 0.0;

   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "first_lp_value", firstdualboundroot) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "first_lp_iterations", scip->stat->nrootfirstlpiterations) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "first_lp_speed", firstlpspeed) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "first_lp_time", firstlptime) );
   SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "final_dual_bound", dualboundroot) );
   SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "final_root_iterations", scip->stat->nrootlpiterations) );

   if( scip->stat->rootlpbestestimate != SCIP_INVALID ) /*lint !e777*/
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "root_lp_estimate", SCIPretransformObj(scip, scip->stat->rootlpbestestimate)) );
   }

   return SCIP_OKAY;
}

/** outputs timing statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_Real readingtime;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintTimingStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   readingtime = SCIPgetReadingTime(scip);

   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "Total Time         : %10.2f\n", readingtime);
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  reading          : %10.2f\n", readingtime);
   }
   else
   {
      SCIP_Real totaltime;
      SCIP_Real solvingtime;

      solvingtime  = SCIPclockGetTime(scip->stat->solvingtime);

      if( scip->set->time_reading )
         totaltime = solvingtime;
      else
         totaltime = solvingtime + readingtime;

      SCIPmessageFPrintInfo(scip->messagehdlr, file, "Total Time         : %10.2f\n", totaltime);
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  solving          : %10.2f\n", solvingtime);
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  presolving       : %10.2f (included in solving)\n", SCIPclockGetTime(scip->stat->presolvingtime));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  reading          : %10.2f%s\n", readingtime, scip->set->time_reading ? " (included in solving)" : "");

      if( scip->stat->ncopies > 0 )
      {
	 SCIP_Real copytime;

	 copytime = SCIPclockGetTime(scip->stat->copyclock);

	 SCIPmessageFPrintInfo(scip->messagehdlr, file, "  copying          : %10.2f (%d #copies) (minimal %.2f, maximal %.2f, average %.2f)\n",
	    copytime, scip->stat->ncopies, scip->stat->mincopytime, scip->stat->maxcopytime, copytime / scip->stat->ncopies);
      }
      else
	 SCIPmessageFPrintInfo(scip->messagehdlr, file, "  copying          : %10.2f %s\n", 0.0, "(0 times copied the problem)");
   }
}

/** collects timing statistics in SCIP_DATATREE
 *
 *  The following keys are set:
 *  - "total_time": Total time spent in SCIP.
 *  - "solving_time": Time spent solving the problem.
 *  - "presolving_time": Time spent in presolving.
 *  - "reading_time": Time spent reading the problem.
 *  - "copy_time": Time spent copying the problem (if applicable).
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPcollectTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_Real readingtime;

   readingtime = SCIPgetReadingTime(scip);
   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
   {
      SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "total_time", readingtime) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "reading_time", readingtime) );
   }
   else
   {
      SCIP_Real totaltime;
      SCIP_Real solvingtime;

      solvingtime = SCIPclockGetTime(scip->stat->solvingtime);

      if( scip->set->time_reading )
         totaltime = solvingtime;
      else
         totaltime = solvingtime + readingtime;
      SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "total_time", totaltime) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "solving_time", solvingtime) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "presolving_time", SCIPclockGetTime(scip->stat->presolvingtime)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "reading_time", readingtime) );

      if( scip->stat->ncopies > 0 )
      {
         SCIP_Real copytime;

         copytime = SCIPclockGetTime(scip->stat->copyclock);
         SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "copy_time", copytime) );
      }
   }

   return SCIP_OKAY;
}

/** outputs expression handler statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintExpressionHandlerStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_Bool headerprinted = FALSE;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintExpressionHandlerStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   for( i = 0; i < scip->set->nexprhdlrs; ++i )
   {
      SCIP_EXPRHDLR* exprhdlr = scip->set->exprhdlrs[i];
      assert(exprhdlr != NULL);

      /* skip unused expression handler */
      if( SCIPexprhdlrGetNCreated(exprhdlr) == 0 )
         continue;

      if( !headerprinted )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file,
            "Expressions        : %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
            "#IntEval", "IntEvalTi", "#RevProp", "RevPropTi", "DomReds", "Cutoffs", "#Estimate", "EstimTime", "Branching", "#Simplify", "SimplifyTi", "Simplified");
         headerprinted = TRUE;
      }

      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17s:", SCIPexprhdlrGetName(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNIntevalCalls(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", SCIPexprhdlrGetIntevalTime(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNReversepropCalls(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", SCIPexprhdlrGetReversepropTime(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNDomainReductions(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNCutoffs(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNEstimateCalls(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", SCIPexprhdlrGetEstimateTime(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNBranchings(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNSimplifyCalls(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", SCIPexprhdlrGetSimplifyTime(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10lld", SCIPexprhdlrGetNSimplifications(exprhdlr));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");
   }
}

/** collects expression handler statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectExpressionHandlerStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* exprhdlrstree;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectExpressionHandlerStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* Create a subtree for expression handler statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &exprhdlrstree, "expression_handlers", scip->set->nexprhdlrs) );

   for( i = 0; i < scip->set->nexprhdlrs; ++i )
   {
      SCIP_DATATREE* exprhdlrdata;
      SCIP_EXPRHDLR* exprhdlr = scip->set->exprhdlrs[i];
      assert(exprhdlr != NULL);

      /* Skip unused expression handlers */
      if( SCIPexprhdlrGetNCreated(exprhdlr) == 0 )
         continue;

      SCIP_CALL( SCIPcreateDatatreeInTree( scip, exprhdlrstree, &exprhdlrdata, SCIPexprhdlrGetName(exprhdlr), 13) );

      SCIP_CALL( SCIPinsertDatatreeString(scip, exprhdlrdata, "description", SCIPexprhdlrGetDescription(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "inteval_calls", SCIPexprhdlrGetNIntevalCalls(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, exprhdlrdata, "inteval_time", SCIPexprhdlrGetIntevalTime(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "reverseprop_calls", SCIPexprhdlrGetNReversepropCalls(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, exprhdlrdata, "reverseprop_time", SCIPexprhdlrGetReversepropTime(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "domain_reductions", SCIPexprhdlrGetNDomainReductions(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "cutoffs", SCIPexprhdlrGetNCutoffs(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "estimate_calls", SCIPexprhdlrGetNEstimateCalls(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, exprhdlrdata, "estimate_time", SCIPexprhdlrGetEstimateTime(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "branchings", SCIPexprhdlrGetNBranchings(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "simplify_calls", SCIPexprhdlrGetNSimplifyCalls(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, exprhdlrdata, "simplify_time", SCIPexprhdlrGetSimplifyTime(exprhdlr)) );
      SCIP_CALL( SCIPinsertDatatreeLong(scip, exprhdlrdata, "simplifications", SCIPexprhdlrGetNSimplifications(exprhdlr)) );
   }

   return SCIP_OKAY;
}

/** outputs NLPI statistics
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
void SCIPprintNLPIStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_Bool printedheader = FALSE;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPprintNLPIStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   for( i = 0; i < scip->set->nnlpis; ++i )
   {
      SCIP_Real solvetime;
      SCIP_Real evaltime = 0.0;
      SCIP_Longint niter;
      SCIP_NLPI* nlpi;
      int j;

      nlpi = scip->set->nlpis[i];
      assert(nlpi != NULL);

      /* skip unused NLP solver */
      if( SCIPnlpiGetNProblems(nlpi) == 0 )
         continue;

      if( !printedheader )
      {
         SCIPmessageFPrintInfo(scip->messagehdlr, file,
            "NLP Solvers        : %10s %10s %10s %10s %s%10s %10s"
            " %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s"
            " %10s %10s %10s %10s %10s %10s %10s\n",
            "#Problems", "ProblemTi", "#Solves", "SolveTime",
            scip->set->time_nlpieval ? "  EvalTime%" : "",
            "#Iter", "Time/Iter",
            "#Okay", "#TimeLimit", "#IterLimit", "#LObjLimit", "#Interrupt", "#NumError", "#EvalError", "#OutOfMem", "#LicenseEr", "#OtherTerm",
            "#GlobOpt", "#LocOpt", "#Feasible", "#LocInfeas", "#GlobInfea", "#Unbounded", "#Unknown"
         );
         printedheader = TRUE;
      }

      solvetime = SCIPnlpiGetSolveTime(nlpi);
      if( scip->set->time_nlpieval )
         evaltime = SCIPnlpiGetEvalTime(nlpi);
      niter = SCIPnlpiGetNIterations(nlpi);

      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  %-17s:", SCIPnlpiGetName(nlpi));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10d", SCIPnlpiGetNProblems(nlpi));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", SCIPnlpiGetProblemTime(nlpi));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10d", SCIPnlpiGetNSolves(nlpi));
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", solvetime);
      if( scip->set->time_nlpieval )
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", solvetime > 0.0 ? 100.0 * evaltime / solvetime : 0.0);
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10" SCIP_LONGINT_FORMAT, niter);
      SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10.2f", niter > 0 ? solvetime / niter : 0.0);

      for( j = (int)SCIP_NLPTERMSTAT_OKAY; j <= (int)SCIP_NLPTERMSTAT_OTHER; ++j )
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10d", SCIPnlpiGetNTermStat(nlpi, (SCIP_NLPTERMSTAT)j));

      for( j = (int)SCIP_NLPSOLSTAT_GLOBOPT; j <= (int)SCIP_NLPSOLSTAT_UNKNOWN; ++j )
         SCIPmessageFPrintInfo(scip->messagehdlr, file, " %10d", SCIPnlpiGetNSolStat(nlpi, (SCIP_NLPSOLSTAT)j));

      SCIPmessageFPrintInfo(scip->messagehdlr, file, "\n");
   }
}

/** give name of NLP termination status as string */
static
const char* nlptermstatToString(
   SCIP_NLPTERMSTAT      termstat            /**< NLP termination status */
   )
{
   switch(termstat)
   {
      case SCIP_NLPTERMSTAT_OKAY:         return "okay";
      case SCIP_NLPTERMSTAT_TIMELIMIT:    return "time_limit";
      case SCIP_NLPTERMSTAT_ITERLIMIT:    return "iter_limit";
      case SCIP_NLPTERMSTAT_LOBJLIMIT:    return "lower_obj_limit";
      case SCIP_NLPTERMSTAT_INTERRUPT:    return "interrupt";
      case SCIP_NLPTERMSTAT_NUMERICERROR: return "numeric_error";
      case SCIP_NLPTERMSTAT_EVALERROR:    return "eval_error";
      case SCIP_NLPTERMSTAT_OUTOFMEMORY:  return "out_of_memory";
      case SCIP_NLPTERMSTAT_LICENSEERROR: return "license_error";
      case SCIP_NLPTERMSTAT_OTHER:        return "other";
      default:                            return "unknown";
   }
}

/** give name of NLP solution status as string */
static
const char* nlpsolstatToString(
   SCIP_NLPSOLSTAT       solstat             /**< NLP solution status */
   )
{
   switch(solstat)
   {
      case SCIP_NLPSOLSTAT_GLOBOPT:         return "global_optimum";
      case SCIP_NLPSOLSTAT_LOCOPT:          return "local_optimum";
      case SCIP_NLPSOLSTAT_FEASIBLE:        return "feasible";
      case SCIP_NLPSOLSTAT_LOCINFEASIBLE:   return "locally_infeasible";
      case SCIP_NLPSOLSTAT_GLOBINFEASIBLE:  return "globally_infeasible";
      case SCIP_NLPSOLSTAT_UNBOUNDED:       return "unbounded";
      case SCIP_NLPSOLSTAT_UNKNOWN:         return "unknown";
      default:                              return "invalid";
   }
}

/** collects NLPI statistics in a SCIP_DATATREE object */
SCIP_RETCODE SCIPcollectNLPIStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_DATATREE* nlpistree;
   int i;

   assert(scip != NULL);
   assert(datatree != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectNLPIStatistics", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->set->nnlpis == 0 )
      return SCIP_OKAY;

   /* Create a subtree for NLPI statistics */
   SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &nlpistree, "nlp_solvers", scip->set->nnlpis) );

   for( i = 0; i < scip->set->nnlpis; ++i )
   {
      SCIP_DATATREE* nlpidata;
      SCIP_Real solvetime;
      SCIP_Real evaltime;
      SCIP_Longint niter;

      SCIP_NLPI* nlpi = scip->set->nlpis[i];
      assert(nlpi != NULL);

      /* Skip unused NLP solvers */
      if( SCIPnlpiGetNProblems(nlpi) == 0 )
         continue;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, nlpistree, &nlpidata, SCIPnlpiGetName(nlpi), 8 + (int)SCIP_NLPTERMSTAT_OTHER + 1 + (int)SCIP_NLPSOLSTAT_UNKNOWN + 1) );

      solvetime = SCIPnlpiGetSolveTime(nlpi);
      evaltime = scip->set->time_nlpieval ? SCIPnlpiGetEvalTime(nlpi) : 0.0;
      niter = SCIPnlpiGetNIterations(nlpi);

      SCIP_CALL( SCIPinsertDatatreeString(scip, nlpidata, "description", SCIPnlpiGetDesc(nlpi)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, nlpidata, "problems", SCIPnlpiGetNProblems(nlpi)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, nlpidata, "problem_time", SCIPnlpiGetProblemTime(nlpi)) );
      SCIP_CALL( SCIPinsertDatatreeInt(scip, nlpidata, "solves", SCIPnlpiGetNSolves(nlpi)) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, nlpidata, "solve_time", solvetime) );
      if( scip->set->time_nlpieval )
      {
         SCIP_CALL( SCIPinsertDatatreeReal(scip, nlpidata, "eval_time_percentage", solvetime > 0.0 ? 100.0 * evaltime / solvetime : 0.0) );
      }
      SCIP_CALL( SCIPinsertDatatreeLong(scip, nlpidata, "iterations", niter) );
      SCIP_CALL( SCIPinsertDatatreeReal(scip, nlpidata, "time_per_iteration", niter > 0 ? solvetime / niter : 0.0) );

      for( int j = (int)SCIP_NLPTERMSTAT_OKAY; j <= (int)SCIP_NLPTERMSTAT_OTHER; ++j )
      {
         SCIP_CALL( SCIPinsertDatatreeInt(scip, nlpidata, nlptermstatToString((SCIP_NLPTERMSTAT)j), SCIPnlpiGetNTermStat(nlpi, (SCIP_NLPTERMSTAT)j)) );
      }

      for( int j = (int)SCIP_NLPSOLSTAT_GLOBOPT; j <= (int)SCIP_NLPSOLSTAT_UNKNOWN; ++j )
      {
         SCIP_CALL( SCIPinsertDatatreeInt(scip, nlpidata, nlpsolstatToString((SCIP_NLPSOLSTAT)j), SCIPnlpiGetNSolStat(nlpi, (SCIP_NLPSOLSTAT)j)) );
      }
   }

   return SCIP_OKAY;
}

/** comparison method for statistics tables */
static
SCIP_DECL_SORTPTRCOMP(tablePosComp)
{  /*lint --e{715}*/
   return (SCIPtableGetPosition((SCIP_TABLE*)elem1) - (SCIPtableGetPosition((SCIP_TABLE*)elem2)));
}

/** outputs solving statistics in JSON format
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note If limits have been changed between the solution and the call to this function, the status is recomputed and
 *        thus may correspond to the original status.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPprintStatisticsJson(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_DATATREE* datatree;
   SCIP_DATATREE* tabledatatree;
   SCIP_TABLE** tables;
   int ntables;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintStatisticsJson", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   ntables = SCIPgetNTables(scip);
   tables = SCIPgetTables(scip);

   /* sort all tables by position unless this has already been done */
   if( !scip->set->tablessorted )
   {
      SCIPsortPtr((void**)tables, tablePosComp, ntables);
      scip->set->tablessorted = TRUE;
   }

   SCIP_CALL( SCIPcreateDatatree(scip, &datatree, ntables) );

   for( i = 0; i < ntables; ++i )
   {
      /* Skip inactive tables or those not relevant to the current stage */
      if( !SCIPtableIsActive(tables[i]) || (SCIPtableGetEarliestStage(tables[i]) > SCIPgetStage(scip)) )
         continue;

      SCIP_CALL( SCIPcreateDatatreeInTree(scip, datatree, &tabledatatree, SCIPtableGetName(tables[i]), -1) );
      SCIP_CALL( SCIPtableCollect(tables[i], scip->set, tabledatatree) );
   }

   SCIP_CALL( SCIPwriteDatatreeJson(scip, file, datatree) );

   SCIPfreeDatatree(scip, &datatree);

   return SCIP_OKAY;
}

/** outputs solving statistics
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @note If limits have been changed between the solution and the call to this function, the status is recomputed and
 *        thus may to correspond to the original status.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPprintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_TABLE** tables;
   int ntables;
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   ntables = SCIPgetNTables(scip);
   tables = SCIPgetTables(scip);

   /* sort all tables by position unless this has already been done */
   if( ! scip->set->tablessorted )
   {
      SCIPsortPtr((void**)tables, tablePosComp, ntables);

      scip->set->tablessorted = TRUE;
   }

   for( i = 0; i < ntables; ++i )
   {
      /* skip tables which are not active or only used in later stages */
      if( ( ! SCIPtableIsActive(tables[i]) ) || SCIPtableGetEarliestStage(tables[i]) > SCIPgetStage(scip) )
         continue;

      SCIP_CALL( SCIPtableOutput(tables[i], scip->mem->probmem, scip->set, file) );
   }

   return SCIP_OKAY;
}

/** outputs reoptimization statistics
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPprintReoptStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_Real solving;
   SCIP_Real presolving;
   SCIP_Real updatetime;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintReoptStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert(scip != NULL);

   /* skip if reoptimization is disabled */
   if( !scip->set->reopt_enable )
      return SCIP_OKAY;

   /* skip if not problem yet */
   if( scip->stat == NULL )
      return SCIP_OKAY;

   solving = SCIPclockGetTime(scip->stat->solvingtimeoverall);
   presolving = SCIPclockGetTime(scip->stat->presolvingtimeoverall);
   updatetime = SCIPclockGetTime(scip->stat->reoptupdatetime);

   SCIPmessageFPrintInfo(scip->messagehdlr, file, "SCIP Reopt Status  : finished after %d runs.\n", scip->stat->nreoptruns);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Time         (sec) :\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  solving          : %10.2f\n", solving);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  presolving       : %10.2f (included in solving)\n", presolving);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  save time        : %10.2f\n", SCIPreoptGetSavingtime(scip->reopt));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  update time      : %10.2f\n", updatetime);
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Nodes              :       feas     infeas     pruned     cutoff\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  total            : %10d %10d %10d %10d\n",
         SCIPreoptGetNTotalFeasNodes(scip->reopt), SCIPreoptGetNTotalInfNodes(scip->reopt),
         SCIPreoptGetNTotalPrunedNodes(scip->reopt), SCIPreoptGetNTotalCutoffReoptnodes(scip->reopt));
   if( scip->stat->nreoptruns > 0 )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  avg              : %10.2f %10.2f %10.2f %10.2f\n",
         (SCIP_Real)SCIPreoptGetNTotalFeasNodes(scip->reopt)/scip->stat->nreoptruns,
         (SCIP_Real)SCIPreoptGetNTotalInfNodes(scip->reopt)/scip->stat->nreoptruns,
         (SCIP_Real)SCIPreoptGetNTotalPrunedNodes(scip->reopt)/scip->stat->nreoptruns,
         (SCIP_Real)SCIPreoptGetNTotalCutoffReoptnodes(scip->reopt)/scip->stat->nreoptruns);
   }
   else
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  avg              : %10s %10s %10s %10s\n", "--", "--", "--", "--");
   }
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "Restarts           :     global      local\n");
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  first            : %10d         --\n", SCIPreoptGetFirstRestarts(scip->reopt));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  last             : %10d         --\n", SCIPreoptGetLastRestarts(scip->reopt));
   SCIPmessageFPrintInfo(scip->messagehdlr, file, "  total            : %10d %10d\n", SCIPreoptGetNRestartsGlobal(scip->reopt),
         SCIPreoptGetNTotalRestartsLocal(scip->reopt));
   if( scip->stat->nreoptruns > 0 )
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  avg              :         -- %10.2f\n",
         (SCIP_Real)SCIPreoptGetNTotalRestartsLocal(scip->reopt)/scip->stat->nreoptruns);
   }
   else
   {
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "  avg              :         -- %10s\n", "--");
   }

   return SCIP_OKAY;
}

/** outputs history statistics about branchings on variables
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
SCIP_RETCODE SCIPprintBranchingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** vars;
   int totalnstrongbranchs;
   int v;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintBranchingStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "problem not yet solved. branching statistics not available.\n");
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, scip->transprob->nvars) );
      for( v = 0; v < scip->transprob->nvars; ++v )
      {
         SCIP_VAR* var;
         int i;

         var = scip->transprob->vars[v];
         for( i = v; i > 0 && strcmp(SCIPvarGetName(var), SCIPvarGetName(vars[i-1])) < 0; i-- )
            vars[i] = vars[i-1];
         vars[i] = var;
      }

      SCIPmessageFPrintInfo(scip->messagehdlr, file, "                                      locks              branchings              inferences      cutoffs                     LP gain          pscostcount                gain variance    \n");
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "variable          prio   factor   down     up  depth    down      up    sb     down       up   down     up            down              up    down      up            down              up\n");

      totalnstrongbranchs = 0;
      for( v = 0; v < scip->transprob->nvars; ++v )
      {
         if( SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS) > 0
            || SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS) > 0
            || SCIPgetVarNStrongbranchs(scip, vars[v]) > 0 )
         {
            int nstrongbranchs;

            nstrongbranchs = SCIPgetVarNStrongbranchs(scip, vars[v]);
            totalnstrongbranchs += nstrongbranchs;
            SCIPmessageFPrintInfo(scip->messagehdlr, file, "%-16s %5d %8.1f %6d %6d %6.1f %7" SCIP_LONGINT_FORMAT " %7" SCIP_LONGINT_FORMAT " %5d %8.1f %8.1f %5.1f%% %5.1f%% %15.4f %15.4f %7.1f %7.1f %15.2f %15.2f\n",
               SCIPvarGetName(vars[v]),
               SCIPvarGetBranchPriority(vars[v]),
               SCIPvarGetBranchFactor(vars[v]),
               SCIPvarGetNLocksDownType(vars[v], SCIP_LOCKTYPE_MODEL),
               SCIPvarGetNLocksUpType(vars[v], SCIP_LOCKTYPE_MODEL),
               (SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_DOWNWARDS)
                  + SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_UPWARDS))/2.0 - 1.0,
               SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS),
               nstrongbranchs,
               SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS),
               100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS),
               100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS),
               SCIPvarGetPseudocost(vars[v], scip->stat, -1.0),
               SCIPvarGetPseudocost(vars[v], scip->stat, +1.0),
               SCIPvarGetPseudocostCount(vars[v], SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetPseudocostCount(vars[v], SCIP_BRANCHDIR_UPWARDS),
               SCIPvarGetPseudocostVariance(vars[v], SCIP_BRANCHDIR_DOWNWARDS, FALSE),
               SCIPvarGetPseudocostVariance(vars[v], SCIP_BRANCHDIR_UPWARDS, FALSE));
         }
      }
      SCIPmessageFPrintInfo(scip->messagehdlr, file, "total                                                %7" SCIP_LONGINT_FORMAT " %7" SCIP_LONGINT_FORMAT " %5d %8.1f %8.1f %5.1f%% %5.1f%% %15.4f %15.4f %7.1f %7.1f %15.2f %15.2f\n",
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS),
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS),
         totalnstrongbranchs,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) > 0
         ? SCIPhistoryGetInferenceSum(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) : 0.0,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) > 0
         ? SCIPhistoryGetInferenceSum(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) : 0.0,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) > 0
         ? SCIPhistoryGetCutoffSum(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) : 0.0,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) > 0
         ? SCIPhistoryGetCutoffSum(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) : 0.0,
         SCIPhistoryGetPseudocost(scip->stat->glbhistory, -1.0),
         SCIPhistoryGetPseudocost(scip->stat->glbhistory, +1.0),
         SCIPhistoryGetPseudocostCount(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS),
         SCIPhistoryGetPseudocostCount(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS),
         SCIPhistoryGetPseudocostVariance(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS),
         SCIPhistoryGetPseudocostVariance(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS));

      SCIPfreeBufferArray(scip, &vars);

      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** collects branching statistics about variables in a SCIP_DATATREE
 *
 * This function collects detailed branching statistics for all variables in the SCIP instance and organizes them into
 * a hierarchical structure in the provided `SCIP_DATATREE`. The statistics include locks, branchings, inferences,
 * cutoffs, pseudocosts, and strong branching information.
 *
 * The `datatree` will contain the following keys:
 * - `variables`: A nested table keyed by variable names, containing:
 *   - `name`: Name of the variable.
 *   - `priority`: Branching priority of the variable.
 *   - `factor`: Branching factor of the variable.
 *   - `locks_down`: Number of locks in the down direction.
 *   - `locks_up`: Number of locks in the up direction.
 *   - `avg_depth`: Average branching depth for the variable.
 *   - `branchings_down`: Number of branchings in the down direction.
 *   - `branchings_up`: Number of branchings in the up direction.
 *   - `strong_branchings`: Number of strong branchings performed on the variable.
 *   - `avg_inferences_down`: Average number of inferences per branching in the down direction.
 *   - `avg_inferences_up`: Average number of inferences per branching in the up direction.
 *   - `cutoff_rate_down`: Percentage of branchings in the down direction that led to cutoffs.
 *   - `cutoff_rate_up`: Percentage of branchings in the up direction that led to cutoffs.
 *   - `pseudocost_down`: Pseudocost in the down direction.
 *   - `pseudocost_up`: Pseudocost in the up direction.
 *   - `pseudocost_count_down`: Number of pseudocost updates in the down direction.
 *   - `pseudocost_count_up`: Number of pseudocost updates in the up direction.
 *   - `pseudocost_variance_down`: Variance of pseudocost in the down direction.
 *   - `pseudocost_variance_up`: Variance of pseudocost in the up direction.
 * - `total_branchings_down`: Total number of branchings in the down direction across all variables.
 * - `total_branchings_up`: Total number of branchings in the up direction across all variables.
 * - `total_strong_branchings`: Total number of strong branchings across all variables.
 * - `avg_inferences_down`: Average inferences per branching in the down direction across all variables.
 * - `avg_inferences_up`: Average inferences per branching in the up direction across all variables.
 * - `avg_cutoff_rate_down`: Average cutoff rate for branchings in the down direction across all variables.
 * - `avg_cutoff_rate_up`: Average cutoff rate for branchings in the up direction across all variables.
 * - `status`: If the problem is not solved, a string indicating that statistics are not available.
 *
 * @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *
 * @return \ref SCIP_OKAY if everything worked. Otherwise, a suitable error code is returned.
 */
SCIP_RETCODE SCIPcollectBranchingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree            /**< data tree */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPcollectBranchingStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
      case SCIP_STAGE_INIT:
      case SCIP_STAGE_PROBLEM:
         SCIP_CALL( SCIPinsertDatatreeString(scip, datatree, "status", "problem not yet solved. branching statistics not available.") );
         return SCIP_OKAY;

      case SCIP_STAGE_TRANSFORMED:
      case SCIP_STAGE_INITPRESOLVE:
      case SCIP_STAGE_PRESOLVING:
      case SCIP_STAGE_EXITPRESOLVE:
      case SCIP_STAGE_PRESOLVED:
      case SCIP_STAGE_SOLVING:
      case SCIP_STAGE_SOLVED:
      {
         SCIP_DATATREE* varsdtree;
         SCIP_VAR** vars;
         int totalnstrongbranchs = 0;
         int v;

         SCIP_CALL( SCIPallocBufferArray(scip, &vars, scip->transprob->nvars) );

         /* Sort variables by name */
         for( v = 0; v < scip->transprob->nvars; ++v )
         {
            SCIP_VAR* var = scip->transprob->vars[v];
            int i;

            for( i = v; i > 0 && strcmp(SCIPvarGetName(var), SCIPvarGetName(vars[i-1])) < 0; i-- )
               vars[i] = vars[i-1];
            vars[i] = var;
         }

         SCIP_CALL( SCIPcreateDatatreeInTree( scip, datatree, &varsdtree, "variables", scip->transprob->nvars + 7 ) );

         /* Collect statistics for each variable */
         for( v = 0; v < scip->transprob->nvars; ++v )
         {
            if( SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS) > 0 ||
               SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS) > 0 ||
               SCIPgetVarNStrongbranchs(scip, vars[v]) > 0 )
            {
               SCIP_DATATREE* vardtree;
               int nstrongbranchs = SCIPgetVarNStrongbranchs(scip, vars[v]);

               totalnstrongbranchs += nstrongbranchs;

               SCIP_CALL( SCIPcreateDatatreeInTree( scip, varsdtree, &vardtree, SCIPvarGetName( vars[v] ), 19 ) );

               SCIP_CALL( SCIPinsertDatatreeString(scip, vardtree, "name", SCIPvarGetName(vars[v])) );
               SCIP_CALL( SCIPinsertDatatreeInt(scip, vardtree, "priority", SCIPvarGetBranchPriority(vars[v])) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "factor", SCIPvarGetBranchFactor(vars[v])) );
               SCIP_CALL( SCIPinsertDatatreeInt(scip, vardtree, "locks_down", SCIPvarGetNLocksDownType(vars[v], SCIP_LOCKTYPE_MODEL)) );
               SCIP_CALL( SCIPinsertDatatreeInt(scip, vardtree, "locks_up", SCIPvarGetNLocksUpType(vars[v], SCIP_LOCKTYPE_MODEL)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "avg_depth",
                  (SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_DOWNWARDS) +
                     SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_UPWARDS)) / 2.0 - 1.0) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, vardtree, "branchings_down", SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeLong(scip, vardtree, "branchings_up", SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeInt(scip, vardtree, "strong_branchings", nstrongbranchs) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "avg_inferences_down", SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "avg_inferences_up", SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "cutoff_rate_down", 100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "cutoff_rate_up", 100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "pseudocost_down", SCIPvarGetPseudocost(vars[v], scip->stat, -1.0)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "pseudocost_up", SCIPvarGetPseudocost(vars[v], scip->stat, +1.0)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "pseudocost_count_down", SCIPvarGetPseudocostCount(vars[v], SCIP_BRANCHDIR_DOWNWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "pseudocost_count_up", SCIPvarGetPseudocostCount(vars[v], SCIP_BRANCHDIR_UPWARDS)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "pseudocost_variance_down", SCIPvarGetPseudocostVariance(vars[v], SCIP_BRANCHDIR_DOWNWARDS, FALSE)) );
               SCIP_CALL( SCIPinsertDatatreeReal(scip, vardtree, "pseudocost_variance_up", SCIPvarGetPseudocostVariance(vars[v], SCIP_BRANCHDIR_UPWARDS, FALSE)) );
            }
         }

         /* add total statistics */
         SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "total_branchings_down", SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS)) );
         SCIP_CALL( SCIPinsertDatatreeLong(scip, datatree, "total_branchings_up", SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS)) );
         SCIP_CALL( SCIPinsertDatatreeInt(scip, datatree, "total_strong_branchings", totalnstrongbranchs) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "avg_inferences_down",
            SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) > 0 ?
            SCIPhistoryGetInferenceSum(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) /
            (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) : 0.0) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "avg_inferences_up",
            SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) > 0 ?
            SCIPhistoryGetInferenceSum(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) /
            (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) : 0.0) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "avg_cutoff_rate_down",
            SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) > 0 ?
            SCIPhistoryGetCutoffSum(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) /
            (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) : 0.0) );
         SCIP_CALL( SCIPinsertDatatreeReal(scip, datatree, "avg_cutoff_rate_up",
            SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) > 0 ?
            SCIPhistoryGetCutoffSum(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) /
            (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) : 0.0) );

         SCIPfreeBufferArray(scip, &vars);
         break;
      }

      default:
         SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
         return SCIP_INVALIDCALL;
   } /*lint !e788*/

   return SCIP_OKAY;
}


/** outputs node information display line
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPprintDisplayLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VERBLEVEL        verblevel,          /**< minimal verbosity level to actually display the information line */
   SCIP_Bool             endline             /**< should the line be terminated with a newline symbol? */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPprintDisplayLine", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   if( (SCIP_VERBLEVEL)scip->set->disp_verblevel >= verblevel )
   {
      SCIP_CALL( SCIPdispPrintLine(scip->set, scip->messagehdlr, scip->stat, file, TRUE, endline) );
   }

   return SCIP_OKAY;
}

/** gets total number of implications between variables that are stored in the implication graph
 *
 *  @return the total number of implications between variables that are stored in the implication graph
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 */
int SCIPgetNImplications(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPgetNImplications", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nimplications;
}

/** stores conflict graph of binary variables' implications into a file, which can be used as input for the DOT tool
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if SCIP is in one of the following stages:
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *
 *  @deprecated because binary implications are now stored as cliques, please use SCIPwriteCliqueGraph() instead
 */ /*lint -e715*/
SCIP_RETCODE SCIPwriteImplicationConflictGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name, or NULL for stdout */
   )
{  /*lint --e{715}*/
   SCIPwarningMessage(scip, "SCIPwriteImplicationConflictGraph() is deprecated and does not do anything anymore. All binary to binary implications are now stored in the clique data structure, which can be written to a GML formatted file via SCIPwriteCliqueGraph().\n");

   return SCIP_OKAY;
}

/** update statistical information when a new solution was found */
void SCIPstoreSolutionGap(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   scip->stat->lastsolgap = SCIPcomputeGap(SCIPsetEpsilon(scip->set), SCIPsetInfinity(scip->set), SCIPgetPrimalbound(scip), SCIPgetDualbound(scip));

   if( scip->primal->nsols == 1 )
      scip->stat->firstsolgap = scip->stat->lastsolgap;

   if( scip->set->misc_calcintegral )
   {
      SCIP_Real upperbound = SCIPgetUpperbound(scip);

      if( upperbound < scip->stat->lastupperbound )
         SCIPstatUpdatePrimalDualIntegrals(scip->stat, scip->set, scip->transprob, scip->origprob, upperbound, -SCIPinfinity(scip));
   }
}

/** recomputes and returns the primal dual gap stored in the stats
 *
 * @return returns the primal dual gap stored in the stats
 */
SCIP_EXPORT
SCIP_Real SCIPgetPrimalDualIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIPstatGetPrimalDualIntegral(scip->stat, scip->set, scip->transprob, scip->origprob, TRUE);
}
