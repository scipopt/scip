/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: var.h,v 1.46 2003/12/15 17:45:35 bzfpfend Exp $"

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
#include "type_event.h"
#include "type_lp.h"
#include "type_var.h"
#include "type_prob.h"
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
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** undoes single bound change */
extern
RETCODE SCIPboundchgUndo(
   BOUNDCHG*        boundchg,           /**< bound change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** frees domain change data */
extern
RETCODE SCIPdomchgFree(
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** converts a dynamic domain change data into a static one, using less memory than for a dynamic one */
extern
RETCODE SCIPdomchgMakeStatic(
   DOMCHG**         domchg,             /**< pointer to domain change data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** applies domain change */
extern
RETCODE SCIPdomchgApply(
   DOMCHG*          domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** undoes domain change */
extern
RETCODE SCIPdomchgUndo(
   DOMCHG*          domchg,             /**< domain change to remove */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** adds bound change to domain changes */
extern
RETCODE SCIPdomchgAddBoundchg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the bounds for */
   Real             newbound,           /**< new value for bound */
   Real             oldbound,           /**< old value for bound */
   BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   NODE*            node,               /**< node where this bound change appears */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   VAR*             infervar            /**< variable that was changed (parent of var, or var itself) */
   );

/** adds hole change to domain changes */
extern
RETCODE SCIPdomchgAddHolechg(
   DOMCHG**         domchg,             /**< pointer to domain change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
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
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
   );

/** creates and captures a loose variable belonging to the transformed problem */
extern
RETCODE SCIPvarCreateTransformed(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             initial,            /**< should var's column be present in the initial root LP? */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
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
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (may be NULL, if it's not a column variable) */
   );

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
extern
RETCODE SCIPvarTransform(
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   OBJSENSE         objsense,           /**< objective sense of original problem; transformed is always MINIMIZE */
   VAR**            transvar            /**< pointer to store the transformed variable */
   );

/** gets corresponding transformed variable of an original or negated original variable */
extern
RETCODE SCIPvarGetTransformed(
   VAR*             origvar,            /**< original problem variable */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR**            transvar            /**< pointer to store the transformed variable, or NULL if not existing yet */
   );

/** converts transformed variable into column variable and creates LP column */
extern
RETCODE SCIPvarColumn(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< actual LP data */
   );

/** converts variable into fixed variable */
extern
RETCODE SCIPvarFix(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< actual LP data */
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
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< actual LP data */
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
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< actual LP data */
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
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newobj              /**< new objective value for variable */
   );

/** adds value to objective value of variable */
extern
RETCODE SCIPvarAddObj(
   VAR*             var,                /**< variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< transformed problem after presolve */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             addobj              /**< additional objective value for variable */
   );

/** changes global lower bound of variable */
extern
RETCODE SCIPvarChgLbGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   );

/** changes global upper bound of variable */
extern
RETCODE SCIPvarChgUbGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound            /**< new bound for variable */
   );

/** changes global bound of variable */
extern
RETCODE SCIPvarChgBdGlobal(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype           /**< type of bound: lower or upper bound */
   );

/** changes current local lower bound of variable */
extern
RETCODE SCIPvarChgLbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   int              inferdepth,         /**< depth in the tree, where this bound change took place */
   int              infernum            /**< bound change index for each node representing the order of changes */
   );

/** changes current local upper bound of variable */
extern
RETCODE SCIPvarChgUbLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   int              inferdepth,         /**< depth in the tree, where this bound change took place */
   int              infernum            /**< bound change index for each node representing the order of changes */
   );

/** changes current local bound of variable */
extern
RETCODE SCIPvarChgBdLocal(
   VAR*             var,                /**< problem variable to change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             newbound,           /**< new bound for variable */
   BOUNDTYPE        boundtype,          /**< type of bound: lower or upper bound */
   CONS*            infercons,          /**< constraint that deduced the bound change (binary variables only), or NULL */
   VAR*             infervar,           /**< variable that was changed (parent of var, or var itself) */
   int              inferdepth,         /**< depth in the tree, where this bound change took place */
   int              infernum            /**< bound change index for each node representing the order of changes */
   );

/** adjust lower bound to integral value, if variable is integral */
extern
void SCIPvarAdjustLb(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            lb                  /**< pointer to lower bound to adjust */
   );

/** adjust upper bound to integral value, if variable is integral */
extern
void SCIPvarAdjustUb(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real*            ub                  /**< pointer to upper bound to adjust */
   );

/** changes lower bound of variable in current dive */
extern
RETCODE SCIPvarChgLbDive(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   );

/** changes upper bound of variable in current dive */
extern
RETCODE SCIPvarChgUbDive(
   VAR*             var,                /**< problem variable to change */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   Real             newbound            /**< new bound for variable */
   );

/** adds a hole to the variable's global domain and to its current local domain */
extern
RETCODE SCIPvarAddHoleGlobal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   );

/** adds a hole to the variable's current local domain */
extern
RETCODE SCIPvarAddHoleLocal(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   Real             left,               /**< left bound of open interval in new hole */
   Real             right               /**< right bound of open interval in new hole */
   );

/** sets the branching priority of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher priority leads to a higher probability that this variable is chosen for branching
 */
extern
void SCIPvarChgBranchingPriority(
   VAR*             var,                /**< problem variable */
   const SET*       set,                /**< global SCIP settings */
   Real             branchingpriority   /**< priority of the variable to choose as branching variable */
   );

/** gets lower bound of variable in current dive */
extern
Real SCIPvarGetLbDive(
   VAR*             var,                /**< problem variable */
   const SET*       set                 /**< global SCIP settings */
   );

/** gets upper bound of variable in current dive */
extern
Real SCIPvarGetUbDive(
   VAR*             var,                /**< problem variable */
   const SET*       set                 /**< global SCIP settings */
   );

/** resolves variable to columns and adds them with the coefficient to the row */
extern
RETCODE SCIPvarAddToRow(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< actual LP data */
   ROW*             row,                /**< LP row */
   Real             val                 /**< value of coefficient */
   );

/** includes event handler with given data in variable's event filter */
extern
RETCODE SCIPvarCatchEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type to catch */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );

/** deletes event handler with given data from variable's event filter */
extern
RETCODE SCIPvarDropEvent(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   EVENTHDLR*       eventhdlr,          /**< event handler to call for the event processing */
   EVENTDATA*       eventdata           /**< event data to pass to the event handler for the event processing */
   );



/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given variable */
extern
DECL_HASHGETKEY(SCIPhashGetKeyVar);



#endif
