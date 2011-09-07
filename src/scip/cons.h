/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons.h
 * @brief  internal methods for constraints and constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_H__
#define __SCIP_CONS_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_mem.h"
#include "scip/type_misc.h"
#include "scip/type_timing.h"
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

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Constraint handler methods
 */

/** copies the given constraint handler to a new scip */
extern
SCIP_RETCODE SCIPconshdlrCopyInclude(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set,                /**< SCIP_SET of SCIP to copy to */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   );

/** creates a constraint handler */
extern
SCIP_RETCODE SCIPconshdlrCreate(
   SCIP_CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of constraint handler */
   const char*           desc,               /**< description of constraint handler */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   int                   enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int                   checkpriority,      /**< priority of the constraint handler for checking feasibility (and propagation) */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int                   eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   SCIP_Bool             delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_Bool             delaypresol,        /**< should presolving method be delayed, if other presolvers found reductions? */
   SCIP_Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagation method of constraint handlers should be executed */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy)),  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   SCIP_DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   SCIP_DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   SCIP_DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   SCIP_DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   SCIP_DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   SCIP_DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   SCIP_DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   SCIP_DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   SCIP_DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   SCIP_DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   SCIP_DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   SCIP_DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   SCIP_DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   SCIP_DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   SCIP_DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   SCIP_DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   SCIP_DECL_CONSCOPY    ((*conscopy)),      /**< constraint copying method */
   SCIP_DECL_CONSPARSE   ((*consparse)),     /**< constraint parsing method */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

/** calls destructor and frees memory of constraint handler */
extern
SCIP_RETCODE SCIPconshdlrFree(
   SCIP_CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls init method of constraint handler */
extern
SCIP_RETCODE SCIPconshdlrInit(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** calls exit method of constraint handler */
extern
SCIP_RETCODE SCIPconshdlrExit(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** informs constraint handler that the presolving process is being started */
extern
SCIP_RETCODE SCIPconshdlrInitpre(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** informs constraint handler that the presolving is finished */
extern
SCIP_RETCODE SCIPconshdlrExitpre(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** informs constraint handler that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPconshdlrInitsol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** informs constraint handler that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPconshdlrExitsol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_Bool             restart             /**< was this exit solve call triggered by a restart? */
   );

/** calls LP initialization method of constraint handler to separate all initial active constraints */
extern
SCIP_RETCODE SCIPconshdlrInitLP(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** calls separator method of constraint handler to separate LP solution */
extern
SCIP_RETCODE SCIPconshdlrSeparateLP(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute separation method even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls separator method of constraint handler to separate given primal solution */
extern
SCIP_RETCODE SCIPconshdlrSeparateSol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute separation method even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for LP solution for all constraints added after last
 *  conshdlrReset() call
 */
extern
SCIP_RETCODE SCIPconshdlrEnforceLPSol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   SCIP_Bool             solinfeasible,      /**< was the solution already found out to be infeasible? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for pseudo solution for all constraints added after last
 *  conshdlrReset() call
 */
extern
SCIP_RETCODE SCIPconshdlrEnforcePseudoSol(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_Bool             solinfeasible,      /**< was the solution already found out to be infeasible? */
   SCIP_Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   SCIP_Bool             forced,             /**< should enforcement of pseudo solution be forced? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls feasibility check method of constraint handler */
extern
SCIP_RETCODE SCIPconshdlrCheck(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls propagation method of constraint handler */
extern
SCIP_RETCODE SCIPconshdlrPropagate(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node; -1 if preprocessing domain propagation */
   SCIP_Bool             fullpropagation,    /**< should all constraints be propagated (or only new ones)? */
   SCIP_Bool             execdelayed,        /**< execute propagation method even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls presolving method of constraint handler */
extern
SCIP_RETCODE SCIPconshdlrPresolve(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_Bool             execdelayed,        /**< execute presolving method even if it is marked to be delayed */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** locks rounding of variables involved in the given constraint constraint handler that doesn't need constraints */
SCIP_RETCODE SCIPconshdlrLockVars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** unlocks rounding of variables involved in the given constraint constraint handler that doesn't need constraints */
extern
SCIP_RETCODE SCIPconshdlrUnlockVars(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   );




/*
 * Constraint set change methods
 */

/** frees constraint set change data and releases all included constraints */
extern
SCIP_RETCODE SCIPconssetchgFree(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** adds constraint addition to constraint set changes, and captures constraint; activates constraint if the
 *  constraint set change data is currently active
 */
extern
SCIP_RETCODE SCIPconssetchgAddAddedCons(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_CONS*            cons,               /**< added constraint */
   int                   depth,              /**< depth of constraint set change's node */
   SCIP_Bool             focusnode,          /**< does the constraint set change belong to the focus node? */
   SCIP_Bool             active              /**< is the constraint set change currently active? */
   );

/** adds constraint disabling to constraint set changes, and captures constraint */
extern
SCIP_RETCODE SCIPconssetchgAddDisabledCons(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS*            cons                /**< disabled constraint */
   );

/** applies constraint set change */
extern
SCIP_RETCODE SCIPconssetchgApply(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of constraint set change's node */
   SCIP_Bool             focusnode           /**< does the constraint set change belong to the focus node? */
   );

/** undoes constraint set change */
extern
SCIP_RETCODE SCIPconssetchgUndo(
   SCIP_CONSSETCHG*      conssetchg,         /**< constraint set change to undo */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** applies constraint set change to the global problem and deletes the constraint set change data */
extern
SCIP_RETCODE SCIPconssetchgMakeGlobal(
   SCIP_CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob                /**< problem data */
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
SCIP_RETCODE SCIPconsCreate(
   SCIP_CONS**           cons,               /**< pointer to constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           name,               /**< name of constraint */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP_CONSDATA*        consdata,           /**< data for this specific constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool             original,           /**< is constraint belonging to the original problem? */
   SCIP_Bool             deleteconsdata      /**< has the constraint data to be deleted if constraint is freed? */
   );

/** copies source constraint of source SCIP into the target constraint for the target SCIP, using the variable map for
 *  mapping the variables of the source SCIP to the variables of the target SCIP; if the copying process was successful
 *  a constraint is creates and captures;
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
SCIP_RETCODE SCIPconsCopy(
   SCIP_CONS**           cons,               /**< pointer to store the created target constraint */
   SCIP_SET*             set,                /**< global SCIP settings of the target SCIP */
   const char*           name,               /**< name of constraint, or NULL if the name of the source constraint should be used */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint handler for this constraint */
   SCIP_CONS*            sourcecons,         /**< source constraint of the source SCIP */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, must not be NULL! */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            success             /**< pointer to store whether the copying was successful or not */
   );

/** parses constraint information (in cip format) out of a string; if the parsing process was successful a constraint is
 *  created, captured, and inserted into the conss array of its constraint handler.
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, an LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
extern
SCIP_RETCODE SCIPconsParse(
   SCIP_CONS**           cons,               /**< pointer to constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   const char*           str,                /**< name of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   );

/** frees a constraint and removes it from the conss array of its constraint handler */
extern
SCIP_RETCODE SCIPconsFree(
   SCIP_CONS**           cons,               /**< constraint to free */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** increases usage counter of constraint */
extern
void SCIPconsCapture(
   SCIP_CONS*            cons                /**< constraint */
   );

/** decreases usage counter of constraint, and frees memory if necessary */
extern
SCIP_RETCODE SCIPconsRelease(
   SCIP_CONS**           cons,               /**< pointer to constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** outputs constraint information to file stream */
extern
SCIP_RETCODE SCIPconsPrint(
   SCIP_CONS*            cons,               /**< constraint to print */
   SCIP_SET*             set,                /**< global SCIP settings */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was created, or from the problem, if it was a problem constraint
 */
extern
SCIP_RETCODE SCIPconsDelete(
   SCIP_CONS*            cons,               /**< constraint to delete */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob                /**< problem data */
   );

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
extern
SCIP_RETCODE SCIPconsTransform(
   SCIP_CONS*            origcons,           /**< original constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory buffer */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   );

/** sets the initial flag of the given constraint */
extern
SCIP_RETCODE SCIPconsSetInitial(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             initial             /**< new value */
   );

/** sets the separate flag of the given constraint */
extern
SCIP_RETCODE SCIPconsSetSeparated(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             separate            /**< new value */
   );

/** sets the enforce flag of the given constraint */
extern
SCIP_RETCODE SCIPconsSetEnforced(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             enforce             /**< new value */
   );

/** sets the check flag of the given constraint */
extern
SCIP_RETCODE SCIPconsSetChecked(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             check               /**< new value */
   );

/** sets the propagate flag of the given constraint */
extern
SCIP_RETCODE SCIPconsSetPropagated(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             propagate           /**< new value */
   );

/** sets the local flag of the given constraint */
extern
void SCIPconsSetLocal(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             local               /**< new value */
   );

/** sets the modifiable flag of the given constraint */
extern
void SCIPconsSetModifiable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             modifiable          /**< new value */
   );

/** sets the dynamic flag of the given constraint */
extern
void SCIPconsSetDynamic(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             dynamic             /**< new value */
   );

/** sets the removable flag of the given constraint */
extern
void SCIPconsSetRemovable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             removable           /**< new value */
   );

/** sets the stickingatnode flag of the given constraint */
extern
void SCIPconsSetStickingAtNode(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             stickingatnode      /**< new value */
   );

/** gives the constraint a new name; ATTENTION: to old pointer is over written that might
 *  result in a memory leakage */
extern
void SCIPconsSetNamePointer(
   SCIP_CONS*            cons,               /**< constraint */
   const char*           name                /**< new name of constraint */
   );

/** gets associated transformed constraint of an original constraint, or NULL if no associated transformed constraint
 *  exists
 */
extern
SCIP_CONS* SCIPconsGetTransformed(
   SCIP_CONS*            cons                /**< constraint */
   );

/** activates constraint or marks constraint to be activated in next update */
extern
SCIP_RETCODE SCIPconsActivate(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth in the tree where the constraint activation takes place, or -1 for global problem */
   SCIP_Bool             focusnode           /**< does the constraint activation take place at the focus node? */
   );

/** deactivates constraint or marks constraint to be deactivated in next update */
extern
SCIP_RETCODE SCIPconsDeactivate(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** enables constraint's separation, enforcing, and propagation capabilities or marks them to be enabled in next update */
extern
SCIP_RETCODE SCIPconsEnable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** disables constraint's separation, enforcing, and propagation capabilities or marks them to be disabled in next update */
extern
SCIP_RETCODE SCIPconsDisable(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/** enables constraint's separation capabilities or marks them to be enabled in next update */
extern
SCIP_RETCODE SCIPconsEnableSeparation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** disables constraint's separation capabilities or marks them to be disabled in next update */
extern
SCIP_RETCODE SCIPconsDisableSeparation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** enables constraint's propagation capabilities or marks them to be enabled in next update */
extern
SCIP_RETCODE SCIPconsEnablePropagation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** disables constraint's propagation capabilities or marks them to be disabled in next update */
extern
SCIP_RETCODE SCIPconsDisablePropagation(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update
 */
extern
SCIP_RETCODE SCIPconsAddAge(
   SCIP_CONS*            cons,               /**< constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_Real             deltaage            /**< value to add to the constraint's age */
   );

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update
 */
extern
SCIP_RETCODE SCIPconsIncAge(
   SCIP_CONS*            cons,               /**< constraint */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob                /**< problem data */
   );

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 *  if it was obsolete, makes constraint useful again or marks constraint to be made useful again in next update
 */
extern
SCIP_RETCODE SCIPconsResetAge(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** resolves the given conflicting bound, that was deduced by the given constraint, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb() and SCIPaddConflictUb()
 */
extern
SCIP_RETCODE SCIPconsResolvePropagation(
   SCIP_CONS*            cons,               /**< constraint that deduced the assignment */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int                   inferinfo,          /**< user inference information attached to the bound change */
   SCIP_BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables */
extern
SCIP_RETCODE SCIPconsAddLocks(
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   nlockspos,          /**< increase in number of rounding locks for constraint */
   int                   nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   );

/** checks single constraint for feasibility of the given solution */
extern
SCIP_RETCODE SCIPconsCheck(
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
extern
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyCons);

#ifdef __cplusplus
}
#endif

#endif
