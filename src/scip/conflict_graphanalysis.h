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

/**@file   conflict_graphanalysis.h
 * @ingroup OTHER_CFILES
 * @brief  methods and datastructures for conflict analysis
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Marc Pfetsch
 * @author Michael Winkler
 * @author Jakob Witzig
 *
 * This file implements a conflict analysis method like the one used in modern
 * SAT solvers like zchaff. The algorithm works as follows:
 *
 * Given is a set of bound changes that are not allowed being applied simultaneously, because they
 * render the current node infeasible (e.g. because a single constraint is infeasible in the these
 * bounds, or because the LP relaxation is infeasible).  The goal is to deduce a clause on variables
 * -- a conflict clause -- representing the "reason" for this conflict, i.e., the branching decisions
 * or the deductions (applied e.g. in domain propagation) that lead to the conflict. This clause can
 * then be added to the constraint set to help cutting off similar parts of the branch and bound
 * tree, that would lead to the same conflict.  A conflict clause can also be generated, if the
 * conflict was detected by a locally valid constraint. In this case, the resulting conflict clause
 * is also locally valid in the same depth as the conflict detecting constraint. If all involved
 * variables are binary, a linear (set covering) constraint can be generated, otherwise a bound
 * disjunction constraint is generated. Details are given in
 *
 * Tobias Achterberg, Conflict Analysis in Mixed Integer Programming@n
 * Discrete Optimization, 4, 4-20 (2007)
 *
 * See also @ref CONF. Here is an outline of the algorithm:
 *
 * -#  Put all the given bound changes to a priority queue, which is ordered,
 *     such that the bound change that was applied last due to branching or deduction
 *     is at the top of the queue. The variables in the queue are always active
 *     problem variables. Because binary variables are preferred over general integer
 *     variables, integer variables are put on the priority queue prior to the binary
 *     variables. Create an empty conflict set.
 * -#  Remove the top bound change b from the priority queue.
 * -#  Perform the following case distinction:
 *     -#  If the remaining queue is non-empty, and bound change b' (the one that is now
 *         on the top of the queue) was applied at the same depth level as b, and if
 *         b was a deduction with known inference reason, and if the inference constraint's
 *         valid depth is smaller or equal to the conflict detecting constraint's valid
 *         depth:
 *          - Resolve bound change b by asking the constraint that inferred the
 *            bound change to put all the bound changes on the priority queue, that
 *            lead to the deduction of b.
 *            Note that these bound changes have at most the same inference depth
 *            level as b, and were deduced earlier than b.
 *     -#  Otherwise, the bound change b was a branching decision or a deduction with
 *         missing inference reason, or the inference constraint's validity is more local
 *         than the one of the conflict detecting constraint.
 *          - If a the bound changed corresponds to a binary variable, add it or its
 *            negation to the conflict set, depending on which of them is currently fixed to
 *            FALSE (i.e., the conflict set consists of literals that cannot be FALSE
 *            altogether at the same time).
 *          - Otherwise put the bound change into the conflict set.
 *         Note that if the bound change was a branching, all deduced bound changes
 *         remaining in the priority queue have smaller inference depth level than b,
 *         since deductions are always applied after the branching decisions. However,
 *         there is the possibility, that b was a deduction, where the inference
 *         reason was not given or the inference constraint was too local.
 *         With this lack of information, we must treat the deduced bound change like
 *         a branching, and there may exist other deduced bound changes of the same
 *         inference depth level in the priority queue.
 * -#  If priority queue is non-empty, goto step 2.
 * -#  The conflict set represents the conflict clause saying that at least one
 *     of the conflict variables must take a different value. The conflict set is then passed
 *     to the conflict handlers, that may create a corresponding constraint (e.g. a logicor
 *     constraint or bound disjunction constraint) out of these conflict variables and
 *     add it to the problem.
 *
 * If all deduced bound changes come with (global) inference information, depending on
 * the conflict analyzing strategy, the resulting conflict set has the following property:
 *  - 1-FirstUIP: In the depth level where the conflict was found, at most one variable
 *    assigned at that level is member of the conflict set. This conflict variable is the
 *    first unique implication point of its depth level (FUIP).
 *  - All-FirstUIP: For each depth level, at most one variable assigned at that level is
 *    member of the conflict set. This conflict variable is the first unique implication
 *    point of its depth level (FUIP).
 *
 * The user has to do the following to get the conflict analysis running in its
 * current implementation:
 *  - A constraint handler or propagator supporting the conflict analysis must implement
 *    the CONSRESPROP/PROPRESPROP call, that processes a bound change inference b and puts all
 *    the reason bounds leading to the application of b with calls to
 *    SCIPaddConflictBound() on the conflict queue (algorithm step 3.(a)).
 *  - If the current bounds lead to a deduction of a bound change (e.g. in domain
 *    propagation), a constraint handler should call SCIPinferVarLbCons() or
 *    SCIPinferVarUbCons(), thus providing the constraint that inferred the bound change.
 *    A propagator should call SCIPinferVarLbProp() or SCIPinferVarUbProp() instead,
 *    thus providing a pointer to itself.
 *  - If (in the current bounds) an infeasibility is detected, the constraint handler or
 *    propagator should
 *     1. call SCIPinitConflictAnalysis() to initialize the conflict queue,
 *     2. call SCIPaddConflictBound() for each bound that lead to the conflict,
 *     3. call SCIPanalyzeConflictCons() or SCIPanalyzeConflict() to analyze the conflict
 *        and add an appropriate conflict constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONFLICT_GRAPHANALYSIS_H__
#define __SCIP_CONFLICT_GRAPHANALYSIS_H__

#include "scip/def.h"
#include "scip/type_cuts.h"
#include "scip/type_conflict.h"
#include "scip/type_reopt.h"
#include "scip/type_implics.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "lpi/type_lpi.h"
#include "scip/type_branch.h"
#include "scip/type_mem.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_event.h"
#include "scip/type_message.h"
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** creates an empty conflict set */
SCIP_RETCODE SCIPconflictsetCreate(
   SCIP_CONFLICTSET**    conflictset,        /**< pointer to store the conflict set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   );

/** frees a conflict set */
void SCIPconflictsetFree(
   SCIP_CONFLICTSET**    conflictset,        /**< pointer to the conflict set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   );

/** copies the given conflict handler to a new scip */
SCIP_RETCODE SCIPconflicthdlrCopyInclude(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a conflict handler */
SCIP_RETCODE SCIPconflicthdlrCreate(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy)),  /**< copy method of conflict handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   );

/** calls destructor and frees memory of conflict handler */
SCIP_RETCODE SCIPconflicthdlrFree(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls init method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs conflict handler that the branch and bound process is being started */
SCIP_RETCODE SCIPconflicthdlrInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs conflict handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPconflicthdlrExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExec(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node to add conflict constraint to */
   SCIP_NODE*            validnode,          /**< node at which the constraint is valid */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change resembling the conflict set */
   SCIP_Real*            relaxedbds,         /**< array with relaxed bounds which are efficient to create a valid conflict */
   int                   nbdchginfos,        /**< number of bound changes in the conflict set */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             usescutoffbound,    /**< depends the conflict on the cutoff bound? */
   SCIP_Bool             resolved,           /**< was the conflict set already used to create a constraint? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of conflict handler */
void SCIPconflicthdlrSetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the conflict handler */
   );

/** set copy method of conflict handler */
void SCIPconflicthdlrSetCopy(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy))   /**< copy method of the conflict handler */
   );

/** set destructor of conflict handler */
void SCIPconflicthdlrSetFree(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTFREE((*conflictfree))   /**< destructor of conflict handler */
   );

/** set initialization method of conflict handler */
void SCIPconflicthdlrSetInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit))   /**< initialization method conflict handler */
   );

/** set deinitialization method of conflict handler */
void SCIPconflicthdlrSetExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit))   /**< deinitialization method conflict handler */
   );

/** set solving process initialization method of conflict handler */
void SCIPconflicthdlrSetInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol))/**< solving process initialization method of conflict handler */
   );

/** set solving process deinitialization method of conflict handler */
void SCIPconflicthdlrSetExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol))/**< solving process deinitialization method of conflict handler */
   );

/** enables or disables all clocks of \p conflicthdlr, depending on the value of the flag */
void SCIPconflicthdlrEnableOrDisableClocks(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< the conflict handler for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the conflict handler be enabled? */
   );

/** return TRUE if conflict analysis is applicable; In case the function return FALSE there is no need to initialize the
 *  conflict analysis since it will not be applied
 */
SCIP_Bool SCIPconflictApplicable(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** creates a temporary bound change information object that is destroyed after the conflict sets are flushed */
SCIP_RETCODE conflictCreateTmpBdchginfo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< active variable that changed the bounds */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BDCHGINFO**      bdchginfo           /**< pointer to store bound change information */
   );


/*
 * Conflict LP Bound Changes
 */


/** calculates the maximal size of conflict sets to be used */
int conflictCalcMaxsize(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   );

/** undoes bound changes on variables, still leaving the given infeasibility proof valid */
SCIP_RETCODE SCIPundoBdchgsProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            proofcoefs,         /**< coefficients in infeasibility proof */
   SCIP_Real             prooflhs,           /**< left hand side of proof */
   SCIP_Real*            proofact,           /**< current activity of proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   SCIP_LPBDCHGS*        oldlpbdchgs,        /**< old LP bound changes used for reset the LP bound change, or NULL */
   SCIP_LPBDCHGS*        relaxedlpbdchgs,    /**< relaxed LP bound changes used for reset the LP bound change, or NULL */
   SCIP_Bool*            resolve,            /**< pointer to store whether the changed LP should be resolved again, or NULL */
   SCIP_LPI*             lpi                 /**< pointer to LPi to access infinity of LP solver; necessary to set correct values */
   );

/** applies conflict analysis starting with given bound changes, that could not be undone during previous
 *  infeasibility analysis
 */
SCIP_RETCODE SCIPconflictAnalyzeRemainingBdchgs(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   );

/** initializes the propagation conflict analysis by clearing the conflict candidate queue */
SCIP_RETCODE SCIPconflictInit(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             usescutoffbound     /**< depends the conflict on a cutoff bound? */
   );

/** adds variable's bound to conflict candidate queue */
SCIP_RETCODE SCIPconflictAddBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   );

/** adds variable's bound to conflict candidate queue with the additional information of a relaxed bound */
SCIP_RETCODE SCIPconflictAddRelaxedBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Real             relaxedbd           /**< the relaxed bound */
   );

/** checks if the given variable is already part of the current conflict set or queued for resolving with the same or
 *  even stronger bound
 */
SCIP_RETCODE SCIPconflictIsVarUsed(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for which the score should be increased */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index (time stamp of bound change), or NULL for current time */
   SCIP_Bool*            used                /**< pointer to store if the variable is already used */
   );

/** returns the conflict lower bound if the variable is present in the current conflict set; otherwise the global lower
 *  bound
 */
SCIP_Real SCIPconflictGetVarLb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** returns the conflict upper bound if the variable is present in the current conflict set; otherwise the global upper
 *  bound
 */
SCIP_Real SCIPconflictGetVarUb(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_VAR*             var                 /**< problem variable */
   );

/** try to find a subset of changed bounds leading to an infeasible LP
 *
 *  1. call undoBdchgsDualfarkas() or undoBdchgsDualsol()
 *     -> update lb/ubchginfoposs arrays
 *     -> store additional changes in bdchg and curvarlbs/ubs arrays
 *     -> apply additional changes to the LPI
 *  2. (optional) if additional bound changes were undone:
 *     -> resolve LP
 *     -> goto 1.
 *  3. redo all bound changes in the LPI to restore the LPI to its original state
 *  4. analyze conflict
 *     -> put remaining changed bounds (see lb/ubchginfoposs arrays) into starting conflict set
 */
SCIP_RETCODE SCIPrunBoundHeuristic(
   SCIP_CONFLICT*        conflict,           /**< conflict data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_LPI*             lpi,                /**< LPI data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Real*            proofcoefs,         /**< coefficients in the proof constraint */
   SCIP_Real*            prooflhs,           /**< lhs of the proof constraint */
   SCIP_Real*            proofactivity,      /**< maximal activity of the proof constraint */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int*                  iterations,         /**< pointer to store the total number of LP iterations used */
   SCIP_Bool             marklpunsolved,     /**< whether LP should be marked unsolved after analysis (needed for strong branching) */
   SCIP_Bool*            dualproofsuccess,   /**< pointer to store success result of dual proof analysis */
   SCIP_Bool*            valid               /**< pointer to store whether the result is still a valid proof */
   );

#ifdef __cplusplus
}
#endif

#endif
