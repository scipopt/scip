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

/**@file   conflict_dualproofanalysis.h
 * @ingroup OTHER_CFILES
 * @brief  internal methods for dual proof conflict analysis
 * @author Timo Berthold
 * @author Jakob Witzig
 *
 * In dual proof analysis, an infeasible LP relaxation is analysed.
 * Using the dual solution, a valid constraint is derived that is violated
 * by all values in the domain. This constraint is added to the problem
 * and can then be used for domain propagation. More details can be found in [1]
 *
 * [1] J. Witzig, T. Berthold, en S. Heinz, ‘Computational aspects of infeasibility analysis in mixed integer programming’,
 * Math. Prog. Comp., mrt. 2021, doi: 10.1007/s12532-021-00202-0.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONFLICT_DUALPROOFANALYSIS_H__
#define __SCIP_CONFLICT_DUALPROOFANALYSIS_H__

#include "blockmemshell/memory.h"
#include "scip/def.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_conflictstore.h"
#include "scip/type_event.h"
#include "scip/type_implics.h"
#include "scip/type_lp.h"
#include "scip/type_prob.h"
#include "scip/type_reopt.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"
#include "scip/type_var.h"
#include "scip/type_cuts.h"

#ifdef __cplusplus
extern "C" {
#endif
/*
 * Proof Sets
 */

/** frees a proofset */
void SCIPproofsetFree(
   SCIP_PROOFSET**       proofset,           /**< proof set */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns the number of variables in the proofset */
int SCIPproofsetGetNVars(
   SCIP_PROOFSET*        proofset            /**< proof set */
   );
/** returns the number of variables in the proofset */


/** creates and clears the proofset */
SCIP_RETCODE SCIPconflictInitProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   );


/* create proof constraints out of proof sets */
SCIP_RETCODE SCIPconflictFlushProofset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict store */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_CLIQUETABLE*     cliquetable         /**< clique table data structure */
   );

/** perform conflict analysis based on a dual unbounded ray
 *
 *  given an aggregation of rows lhs <= a^Tx such that lhs > maxactivity. if the constraint has size one we add a
 *  bound change instead of the constraint.
 */
SCIP_RETCODE SCIPconflictAnalyzeDualProof(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< tree data */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_AGGRROW*         proofrow,           /**< aggregated row representing the proof */
   int                   validdepth,         /**< valid depth of the dual proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   SCIP_Bool             initialproof,       /**< do we analyze the initial reason of infeasibility? */
   SCIP_Bool*            globalinfeasible,   /**< pointer to store whether global infeasibility could be proven */
   SCIP_Bool*            success             /**< pointer to store success result */
   );

#ifdef __cplusplus
}
#endif

#endif
