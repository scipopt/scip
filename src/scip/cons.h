/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons.h,v 1.90 2005/02/28 13:26:22 bzfpfend Exp $"

/**@file   cons.h
 * @brief  internal methods for constraints and constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_H__
#define __CONS_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_mem.h"
#include "scip/type_misc.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_tree.h"
#include "scip/type_sepastore.h"
#include "scip/type_cons.h"
#include "scip/pub_cons.h"

#ifndef NDEBUG
#include "scip/struct_cons.h"
#endif



/*
 * Constraint handler methods
 */

/** creates a constraint handler */
extern
RETCODE SCIPconshdlrCreate(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              checkpriority,      /**< priority of the constraint handler for checking feasibility */
   int              sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int              eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int              maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   Bool             delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
   Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   Bool             delaypresol,        /**< should presolving method be delayed, if other presolvers found reductions? */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   DECL_CONSSEPA    ((*conssepa)),      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

/** calls destructor and frees memory of constraint handler */
extern
RETCODE SCIPconshdlrFree(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SET*             set                 /**< global SCIP settings */
   );

/** calls init method of constraint handler */
extern
RETCODE SCIPconshdlrInit(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** calls exit method of constraint handler */
extern
RETCODE SCIPconshdlrExit(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** informs constraint handler that the presolving process is being started */
extern
RETCODE SCIPconshdlrInitpre(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** informs constraint handler that the presolving is finished */
extern
RETCODE SCIPconshdlrExitpre(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** informs constraint handler that the branch and bound process is being started */
extern
RETCODE SCIPconshdlrInitsol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** informs constraint handler that the branch and bound process data is being freed */
extern
RETCODE SCIPconshdlrExitsol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** calls LP initialization method of constraint handler to separate all initial active constraints */
extern
RETCODE SCIPconshdlrInitLP(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** calls separator method of constraint handler to separate all constraints added after last conshdlrReset() call */
extern
RETCODE SCIPconshdlrSeparate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              depth,              /**< depth of current node */
   Bool             execdelayed,        /**< execute separation method even if it is marked to be delayed */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for LP solution for all constraints added after last
 *  conshdlrReset() call
 */
extern
RETCODE SCIPconshdlrEnforceLPSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   SEPASTORE*       sepastore,          /**< separation storage */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for pseudo solution for all constraints added after last
 *  conshdlrReset() call
 */
extern
RETCODE SCIPconshdlrEnforcePseudoSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls feasibility check method of constraint handler */
extern
RETCODE SCIPconshdlrCheck(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls propagation method of constraint handler */
extern
RETCODE SCIPconshdlrPropagate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth,              /**< depth of current node; -1 if preprocessing domain propagation */
   Bool             execdelayed,        /**< execute propagation method even if it is marked to be delayed */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls presolving method of constraint handler */
extern
RETCODE SCIPconshdlrPresolve(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   Bool             execdelayed,        /**< execute presolving method even if it is marked to be delayed */
   int              nrounds,            /**< number of presolving rounds already done */
   int*             nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*             naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*             nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*             nchgbds,            /**< pointer to total number of variable bounds tightend of all presolvers */
   int*             naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*             ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*             nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*             nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*             nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** locks rounding of variables involved in the given constraint constraint handler that doesn't need constraints */
RETCODE SCIPconshdlrLockVars(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set                 /**< global SCIP settings */
   );

/** unlocks rounding of variables involved in the given constraint constraint handler that doesn't need constraints */
extern
RETCODE SCIPconshdlrUnlockVars(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   SET*             set                 /**< global SCIP settings */
   );




/*
 * Constraint set change methods
 */

/** frees constraint set change data and releases all included constraints */
extern
RETCODE SCIPconssetchgFree(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** adds constraint addition to constraint set changes, and captures constraint; activates constraint if the
 *  constraint set change data is currently active
 */
extern
RETCODE SCIPconssetchgAddAddedCons(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   CONS*            cons,               /**< added constraint */
   int              depth,              /**< depth of constraint set change's node */
   Bool             active              /**< is the constraint set change currently active? */
   );

/** adds constraint disabling to constraint set changes, and captures constraint */
extern
RETCODE SCIPconssetchgAddDisabledCons(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   CONS*            cons                /**< disabled constraint */
   );

/** applies constraint set change */
extern
RETCODE SCIPconssetchgApply(
   CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth               /**< depth of constraint set change's node */
   );

/** undoes constraint set change */
extern
RETCODE SCIPconssetchgUndo(
   CONSSETCHG*      conssetchg,         /**< constraint set change to undo */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );




/*
 * Constraint methods
 */

/** creates and captures a constraint, and inserts it into the conss array of its constraint handler
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
extern
RETCODE SCIPconsCreate(
   CONS**           cons,               /**< pointer to constraint */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable,         /**< should the constraint be removed from the LP due to aging or cleanup? */
   Bool             original            /**< is constraint belonging to the original problem? */
   );

/** frees a constraint and removes it from the conss array of its constraint handler */
extern
RETCODE SCIPconsFree(
   CONS**           cons,               /**< constraint to free */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set                 /**< global SCIP settings */
   );

/** increases usage counter of constraint */
extern
void SCIPconsCapture(
   CONS*            cons                /**< constraint */
   );

/** decreases usage counter of constraint, and frees memory if necessary */
extern
RETCODE SCIPconsRelease(
   CONS**           cons,               /**< pointer to constraint */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** outputs constraint information to file stream */
extern
RETCODE SCIPconsPrint(
   CONS*            cons,               /**< constraint to print */
   SET*             set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was created, or from the problem, if it was a problem constraint
 */
extern
RETCODE SCIPconsDelete(
   CONS*            cons,               /**< constraint to delete */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob                /**< problem data */
   );

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
extern
RETCODE SCIPconsTransform(
   CONS*            origcons,           /**< original constraint */
   BLKMEM*          blkmem,             /**< block memory buffer */
   SET*             set,                /**< global SCIP settings */
   CONS**           transcons           /**< pointer to store the transformed constraint */
   );

/** marks the constraint to be globally valid */
extern
void SCIPconsSetGlobal(
   CONS*            cons                /**< constraint */
   );

/** gets associated transformed constraint of an original constraint, or NULL if no associated transformed constraint
 *  exists
 */
extern
CONS* SCIPconsGetTransformed(
   CONS*            cons                /**< constraint */
   );

/** activates constraint or marks constraint to be activated in next update */
extern
RETCODE SCIPconsActivate(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth               /**< depth in the tree where the constraint activation takes place, or -1 for global problem */
   );

/** deactivates constraint or marks constraint to be deactivated in next update */
extern
RETCODE SCIPconsDeactivate(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** enables constraint's separation, enforcing, and propagation capabilities or marks them to be enabled in next update */
extern
RETCODE SCIPconsEnable(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** disables constraint's separation, enforcing, and propagation capabilities or marks them to be disabled in next update */
extern
RETCODE SCIPconsDisable(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat                /**< dynamic problem statistics */
   );

/** enables constraint's propagation capabilities or marks them to be enabled in next update */
extern
RETCODE SCIPconsEnablePropagation(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   );

/** disables constraint's propagation capabilities or marks them to be disabled in next update */
extern
RETCODE SCIPconsDisablePropagation(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   );

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update,
 */
extern
RETCODE SCIPconsAddAge(
   CONS*            cons,               /**< constraint */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< problem data */
   Real             deltaage            /**< value to add to the constraint's age */
   );

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update,
 */
extern
RETCODE SCIPconsIncAge(
   CONS*            cons,               /**< constraint */
   BLKMEM*          blkmem,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob                /**< problem data */
   );

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 *  if it was obsolete, makes constraint useful again or marks constraint to be made useful again in next update
 */
extern
RETCODE SCIPconsResetAge(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   );

/** resolves the given conflicting bound, that was deduced by the given constraint, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb() and SCIPaddConflictUb()
 */
extern
RETCODE SCIPconsResolvePropagation(
   CONS*            cons,               /**< constraint that deduced the assignment */
   SET*             set,                /**< global SCIP settings */
   VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int              inferinfo,          /**< user inference information attached to the bound change */
   BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables */
extern
RETCODE SCIPconsAddLocks(
   CONS*            cons,               /**< constraint */
   SET*             set,                /**< global SCIP settings */
   int              nlockspos,          /**< increase in number of rounding locks for constraint */
   int              nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   );

/** checks single constraint for feasibility of the given solution */
extern
RETCODE SCIPconsCheck(
   CONS*            cons,               /**< constraint to check */
   SET*             set,                /**< global SCIP settings */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** marks the constraint to be essential for feasibility */
extern
RETCODE SCIPconsSetChecked(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   );




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
extern
DECL_HASHGETKEY(SCIPhashGetKeyCons);


#endif
