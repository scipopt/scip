/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: var.h,v 1.65 2004/05/03 13:35:25 bzfpfend Exp $"

/**@file   var.h
 * @brief  internal methods for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __VAR_H__
#define __VAR_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_misc.h"
#include "type_history.h"
#include "type_event.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
#include "type_primal.h"
#include "type_branch.h"
#include "type_cons.h"
#include "pub_var.h"

#ifndef NDEBUG
#include "struct_var.h"
#endif




/*
 * domain change methods
 */

/** applies single bound change */
extern
RETCODE SCIPboundchgApply(
   BOUNDCHG*        boundchg,           /**< bound change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              depth,              /**< depth in the tree, where the bound change takes place */
   int              index               /**< index of the bound change in the bound change array */
   );

/** undoes single bound change */
extern
RETCODE SCIPboundchgUndo(
   BOUNDCHG*        boundchg,           /**< bound change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** frees domain change data */
extern
RETCODE SCIPdomchgFree(
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** converts a dynamic domain change data into a static one, using less memory than for a dynamic one */
extern
RETCODE SCIPdomchgMakeStatic(
   DOMCHG**         domchg,             /**< pointer to domain change data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set                 /**< global SCIP settings */
   );

/** applies domain change */
extern
RETCODE SCIPdomchgApply(
   DOMCHG*          domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              depth               /**< depth in the tree, where the domain change takes place */
   );

/** undoes domain change */
extern
RETCODE SCIPdomchgUndo(
   DOMCHG*          domchg,             /**< domain change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** adds bound change to domain changes */
extern
RETCODE SCIPdomchgAddBoundchg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   BOUNDCHGTYPE     boundchgtype,       /**< type of bound change: branching decision or inference */
   Real             lpsolval,           /**< solval of variable in last LP on path to node, or SCIP_INVALID if unknown */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo           /**< user information for inference to help resolving the conflict */
   );

/** adds hole change to domain changes */
extern
RETCODE SCIPdomchgAddHolechg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   HOLELIST**       ptr,                /**< changed list pointer */
   HOLELIST*        newlist,            /**< new value of list pointer */
   HOLELIST*        oldlist             /**< old value of list pointer */
   );




/*
 * methods for variables 
 */

/** creates and captures an original problem variable */
extern
RETCODE SCIPvarCreateOriginal(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable,         /**< is var's column removeable from the LP (due to aging or cleanup)? */
   DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   VARDATA*         vardata             /**< user data for this specific variable */
   );

/** creates and captures a loose variable belonging to the transformed problem */
extern
RETCODE SCIPvarCreateTransformed(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable,         /**< is var's column removeable from the LP (due to aging or cleanup)? */
   DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   VARDATA*         vardata             /**< user data for this specific variable */
   );

/** increases usage counter of variable */
extern
void SCIPvarCapture(
   VAR*             var                 /**< variable */
   );

/** decreases usage counter of variable, and frees memory if necessary */
extern
RETCODE SCIPvarRelease(
   VAR**            var,                /**< pointer to variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   LP*              lp                  /**< current LP data (may be NULL, if it's not a column variable) */
   );

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
extern
RETCODE SCIPvarTransform(
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets corresponding transformed variable of an original or negated original variable */
extern
RETCODE SCIPvarGetTransformed(
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR**            transvar            /**< pointer to store the transformed variable, or NULL if not existing yet */
   );

/** converts transformed variable into column variable and creates LP column */
extern
RETCODE SCIPvarColumn(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< current LP data */
   );

/** converts column transformed variable back into loose variable, frees LP column */
extern
RETCODE SCIPvarLoose(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< current LP data */
   );

/** converts variable into fixed variable */
extern
RETCODE SCIPvarFix(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             fixedval,           /**< value to fix variable at */
   Bool*            infeasible          /**< pointer to store whether the fixing is infeasible */
   );

/** converts loose variable into aggregated variable */
extern
RETCODE SCIPvarAggregate(
   VAR*             var,                /**< loose problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             aggvar,             /**< loose variable y in aggregation x = a*y + c */
   Real             scalar,             /**< multiplier a in aggregation x = a*y + c */
   Real             constant,           /**< constant shift c in aggregation x = a*y + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   );

/** converts variable into multi-aggregated variable */
extern
RETCODE SCIPvarMultiaggregate(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   );

/** gets negated variable x' = offset - x of problem variable x; the negated variable is created if not yet existing;
 *  the negation offset of binary variables is always 1, the offset of other variables is fixed to lb + ub when the
 *  negated variable is created
 */
extern
RETCODE SCIPvarNegate(
   VAR*             var,                /**< problem variable to negate */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR**            negvar              /**< pointer to store the negated variable */
   );

/** changes type of variable; cannot be called, if var belongs to a problem */
extern
RETCODE SCIPvarChgType(
   VAR*             var,                /**< problem variable to change */
   VARTYPE          vartype             /**< new type of variable */
   );

/** changes objective value of variable */
extern
RETCODE SCIPvarChgObj(
   VAR*             var,                /**< variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newobj              /**< new objective value for variable */
   );

/** adds value to objective value of variable */
extern
RETCODE SCIPvarAddObj(
   VAR*             var,                /**< variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   PRIMAL*          primal,             /**< primal data */
   LP*              lp,                 /**< current LP data */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             addobj              /**< additional objective value for variable */
   );

/** changes objective value of variable in current dive */
extern
RETCODE SCIPvarChgObjDive(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newobj              /**< new objective value for variable */
   );

/** adjust lower bound to integral value, if variable is integral */
extern
void SCIPvarAdjustLb(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real*            lb                  /**< pointer to lower bound to adjust */
   );

/** adjust upper bound to integral value, if variable is integral */
extern
void SCIPvarAdjustUb(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real*            ub                  /**< pointer to upper bound to adjust */
   );

/** sets global lower bound of variable; if possible, adjusts bound to integral value */
extern
RETCODE SCIPvarSetLbGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   );

/** sets global upper bound of variable; if possible, adjusts bound to integral value */
extern
RETCODE SCIPvarSetUbGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   );

/** sets global bound of variable; if possible, adjusts bound to integral value */
extern
RETCODE SCIPvarSetBdGlobal(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** changes current local lower bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
extern
RETCODE SCIPvarChgLbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place, or -1 */
   int              fixindex,           /**< bound change index for each node representing the order of changes, or -1 */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   );

/** changes current local upper bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
extern
RETCODE SCIPvarChgUbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place, or -1 */
   int              fixindex,           /**< bound change index for each node representing the order of changes, or -1 */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   );

/** changes current local bound of variable; if possible, adjusts bound to integral value; stores inference
 *  information in variable
 */
extern
RETCODE SCIPvarChgBdLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself), or NULL */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   int              inferinfo,          /**< user information for inference to help resolving the conflict */
   int              fixdepth,           /**< depth in the tree, where this bound change took place, or -1 */
   int              fixindex,           /**< bound change index for each node representing the order of changes, or -1 */
   BOUNDCHGTYPE     boundchgtype        /**< bound change type (branching or inference) of binary variable's fixing */
   );

/** sets current local lower bound of variable; if possible, adjusts bound to integral value; doesn't change the
 *  inference information stored in the variable
 */
extern
RETCODE SCIPvarSetLbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   );

/** sets current local upper bound of variable; if possible, adjusts bound to integral value; doesn't change the
 *  inference information stored in the variable
 */
extern
RETCODE SCIPvarSetUbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound            /**< new bound for variable */
   );

/** sets current local bound of variable; if possible, adjusts bound to integral value; doesn't change the
 *  inference information stored in the variable
 */
extern
RETCODE SCIPvarSetBdLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** changes lower bound of variable in current dive; if possible, adjusts bound to integral value */
extern
RETCODE SCIPvarChgLbDive(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newbound            /**< new bound for variable */
   );

/** changes upper bound of variable in current dive; if possible, adjusts bound to integral value */
extern
RETCODE SCIPvarChgUbDive(
   VAR*             var,                /**< problem variable to change */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< current LP data */
   Real             newbound            /**< new bound for variable */
   );

/** adds a hole to the variable's global domain and to its current local domain */
extern
RETCODE SCIPvarAddHoleGlobal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   );

/** adds a hole to the variable's current local domain */
extern
RETCODE SCIPvarAddHoleLocal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   );

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z */
extern
RETCODE SCIPvarAddVlb(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   Real             vlbconstant         /**< constant d    in x >= b*z + d */
   );

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z */
extern
RETCODE SCIPvarAddVub(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   Real             vubconstant         /**< constant d    in x <= b*z + d */
   );

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
extern
void SCIPvarChgBranchFactor(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   Real             branchfactor        /**< factor to weigh variable's branching score with */
   );

/** sets the branch priority of the variable; variables with higher branch priority are always prefered to variables
 *  with lower priority in selection of branching variable
 */
extern
void SCIPvarChgBranchPriority(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   int              branchpriority      /**< branching priority of the variable */
   );

/** gets objective value of variable in current dive */
extern
Real SCIPvarGetObjDive(
   VAR*             var,                /**< problem variable */
   SET*             set                 /**< global SCIP settings */
   );

/** gets lower bound of variable in current dive */
extern
Real SCIPvarGetLbDive(
   VAR*             var,                /**< problem variable */
   SET*             set                 /**< global SCIP settings */
   );

/** gets upper bound of variable in current dive */
extern
Real SCIPvarGetUbDive(
   VAR*             var,                /**< problem variable */
   SET*             set                 /**< global SCIP settings */
   );

/** resolves variable to columns and adds them with the coefficient to the row */
extern
RETCODE SCIPvarAddToRow(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< current LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

/** includes event handler with given data in variable's event filter */
extern
RETCODE SCIPvarCatchEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );

/** deletes event handler with given data from variable's event filter */
extern
RETCODE SCIPvarDropEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of
 *  "solvaldelta" in the variable's solution value and resulting change of "objdelta" in the in the LP's objective value
 */
extern
RETCODE SCIPvarUpdatePseudocost(
   VAR*             var,                /**< problem variable */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   );

/** gets the variable's pseudo cost value for the given step size "solvaldelta" in the variable's LP solution value */
extern
Real SCIPvarGetPseudocost(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   );

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
extern
Real SCIPvarGetPseudocostCount(
   VAR*             var,                /**< problem variable */
   BRANCHDIR        dir                 /**< branching direction: 0 (down), or 1 (up) */
   );

/** increases the number of branchings counter of the variable */
extern
RETCODE SCIPvarIncNBranchings(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   int              depth,              /**< depth at which the bound change took place */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** increases the number of inferences counter of the variable */
extern
RETCODE SCIPvarIncNInferences(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction */
   );

/** returns the average number of inferences found after branching on the variable in given direction */
extern
Real SCIPvarGetAvgInferences(
   VAR*             var,                /**< problem variable */
   STAT*            stat,               /**< problem statistics */
   BRANCHDIR        dir                 /**< branching direction */
   );




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
extern
DECL_HASHGETKEY(SCIPhashGetKeyVar);



#endif
