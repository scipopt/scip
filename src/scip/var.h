/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   var.h
 * @brief  Methods and datastructures for problem variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __VAR_H__
#define __VAR_H__



/** status of problem variables */
enum Varstatus
{
   SCIP_VARSTATUS_ORIGINAL   = 0,       /**< variable belongs to original problem */
   SCIP_VARSTATUS_LOOSE      = 1,       /**< variable is a loose variable of the transformed problem */
   SCIP_VARSTATUS_COLUMN     = 2,       /**< variable is a column of the transformed problem */
   SCIP_VARSTATUS_FIXED      = 3,       /**< variable is fixed to specific value in the transformed problem */
   SCIP_VARSTATUS_AGGREGATED = 4,       /**< variable is aggregated to $x = a*y + c$ in the transformed problem */
   SCIP_VARSTATUS_MULTAGGR   = 5,       /**< variable is aggregated to $x = a_1*y_1 + ... + a_k*y_k + c$ */
   SCIP_VARSTATUS_NEGATED    = 6        /**< variable is the negation of an original or transformed variable */
};
typedef enum Varstatus VARSTATUS;

/** variable type */
enum Vartype
{
   SCIP_VARTYPE_BINARY    = 0,          /**< binary variable: $x \in \{0,1\}$ */
   SCIP_VARTYPE_INTEGER   = 1,          /**< integer variable: $x \in \{lb, \ldots, \ub\}$ */
   SCIP_VARTYPE_IMPLINT   = 2,          /**< implicit integer variable: continous variable, that is always integral */
   SCIP_VARTYPE_CONTINOUS = 3           /**< continous variable: $x \in [lb,ub] */
};
typedef enum Vartype VARTYPE;

/** domain change data type */
enum DomchgType
{
   SCIP_DOMCHGTYPE_DYNAMIC = 0,         /**< dynamic bound changes with size information of arrays */
   SCIP_DOMCHGTYPE_BOTH    = 1,         /**< static domain changes: number of entries equals size of arrays */
   SCIP_DOMCHGTYPE_BOUND   = 2          /**< static domain changes without any hole changes */
};
typedef enum DomchgType DOMCHGTYPE;

typedef struct DomChgBound DOMCHGBOUND; /**< static domain change for bound changes */
typedef struct DomChgBoth DOMCHGBOTH;   /**< static domain change for bound and hole changes */
typedef struct DomChgDyn DOMCHGDYN;     /**< dynamic domain change for bound and hole changes */
typedef union DomChg DOMCHG;            /**< changes in domains of variables */
typedef struct BoundChg BOUNDCHG;       /**< changes in bounds of variables */
typedef struct HoleChg HOLECHG;         /**< changes in holelist of variables */
typedef struct Hole HOLE;               /**< hole in a domain of an integer variable */
typedef struct Holelist HOLELIST;       /**< list of holes in a domain of an integer variable */
typedef struct Dom DOM;                 /**< datastructures for storing domains of variables */
typedef struct Aggregate AGGREGATE;     /**< aggregation information */
typedef struct Multaggr MULTAGGR;       /**< multiple aggregation information */
typedef struct Negate NEGATE;           /**< negation information */
typedef struct Var VAR;                 /**< variable of the problem */



#include "def.h"
#include "memory.h"
#include "retcode.h"
#include "lp.h"
#include "stat.h"
#include "tree.h"
#include "prob.h"
#include "cons.h"
#include "branch.h"
#include "event.h"


/** hole in a domain */
struct Hole
{
   Real             left;               /**< left bound of open interval defining the hole $(left,right)$ */
   Real             right;              /**< right bound of open interval defining the hole $(left,right)$ */
};

/** list of domain holes */
struct Holelist
{
   HOLE             hole;               /**< this hole */
   HOLELIST*        next;               /**< next hole in list */
};

/** change in a hole list */
struct HoleChg
{
   HOLELIST**       ptr;                /**< changed list pointer */
   HOLELIST*        newlist;            /**< new value of list pointer */
   HOLELIST*        oldlist;            /**< old value of list pointer */
};

/** change in one bound of a variable */
struct BoundChg
{
   VAR*             var;                /**< active variable to change the bounds for */
   CONS*            infercons;          /**< constraint that infered this bound change, or NULL */
   VAR*             infervar;           /**< variable that was changed (parent of var, or var itself) */
   Real             newbound;           /**< new value for bound */
   Real             oldbound;           /**< old value for bound */
   unsigned int     inferdepth:16;      /**< depth in the tree, where this bound change took place */
   unsigned int     infernum:15;        /**< bound change index for each node representing the order of changes */
   unsigned int     boundtype:1;        /**< type of bound for var: lower or upper bound */
};

/** tracks changes of the variable's domains (static arrays, bound changes only) */
struct DomChgBound
{
   unsigned int     domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   unsigned int     nboundchgs:30;      /**< number of bound changes */
   BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
};

/** tracks changes of the variable's domains (static arrays, bound and hole changes) */
struct DomChgBoth
{
   unsigned int     domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   unsigned int     nboundchgs:30;      /**< number of bound changes */
   BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
   int              nholechgs;          /**< number of hole list changes */
   HOLECHG*         holechgs;           /**< array with changes in hole lists */
};

/** tracks changes of the variable's domains (dynamic arrays) */
struct DomChgDyn
{
   unsigned int     domchgtype:2;       /**< type of domain change data (must be first structure entry!) */
   unsigned int     nboundchgs:30;      /**< number of bound changes */
   BOUNDCHG*        boundchgs;          /**< array with changes in bounds of variables */
   int              nholechgs;          /**< number of hole list changes */
   HOLECHG*         holechgs;           /**< array with changes in hole lists */
   int              boundchgssize;      /**< size of bound changes array */
   int              holechgssize;       /**< size of hole changes array */
};

/** tracks changes of the variable's domains */
union DomChg
{
   DOMCHGBOUND      domchgbound;        /**< bound changes */
   DOMCHGBOTH       domchgboth;         /**< bound and hole changes */
   DOMCHGDYN        domchgdyn;          /**< bound and hole changes with dynamic arrays */
};

/** domain of a variable */
struct Dom
{
   HOLELIST*        holelist;           /**< list of holes */
   Real             lb;                 /**< lower bounds of variables */
   Real             ub;                 /**< upper bounds of variables */
};

/** aggregation information: $x = a*y + c$ */
struct Aggregate
{
   VAR*             var;                /**< variable $y$ in aggregation */
   Real             scalar;             /**< multiplier $a$ in aggregation */
   Real             constant;           /**< constant shift $c$ in aggregation */
};

/** multiple aggregation information: $x = a_1*y_1 + ... + a_k*y_k + c$ */
struct Multaggr
{
   VAR**            vars;               /**< variables $y$ in multiple aggregation */
   Real*            scalars;            /**< multipliers $a$ in multiple aggregation */
   Real             constant;           /**< constant shift $c$ in multiple aggregation */
   int              nvars;              /**< number of variables in aggregation */
   int              varssize;           /**< size of vars and scalars arrays */
};

/** negation information: $x' = c - x$ */
struct Negate
{
   Real             constant;           /**< constant shift $c$ in negation */
};

/** variable of the problem */
struct Var
{
   union
   {
      VAR*          transvar;           /**< pointer to representing transformed variable (for original variables) */
      COL*          col;                /**< LP column (for column variables) */
      AGGREGATE     aggregate;          /**< aggregation information (for aggregated variables) */
      MULTAGGR      multaggr;           /**< multiple aggregation information (for multiple aggregated variables) */
      NEGATE        negate;             /**< negation information (for negated variables) */
   } data;
   char*            name;               /**< name of the variable */
   VAR**            parentvars;         /**< parent variables in the aggregation tree */
   VAR*             negatedvar;         /**< pointer to the variables negation: x' = lb + ub - x, or NULL if not created */
   EVENTFILTER*     eventfilter;        /**< event filter for events concerning this variable; not for ORIGINAL vars */
   CONS*            infercons;          /**< constraint that deduced the assignment (binary variables only), or NULL */
   VAR*             infervar;           /**< variable that was assigned (parent of var, or var itself) */
   DOM              glbdom;             /**< domain of variable in global problem */
   DOM              actdom;             /**< domain of variable in actual subproblem */
   Real             obj;                /**< objective function value of variable */
   int              index;              /**< consecutively numbered variable identifier */
   int              probindex;          /**< array position in problems vars array, or -1 if not assigned to a problem */
   int              pseudocandindex;    /**< array position in pseudo branching candidates array, or -1 */
   int              eventqueueindexobj; /**< array position in event queue of objective change event, or -1 */
   int              eventqueueindexlb;  /**< array position in event queue of lower bound change event, or -1 */
   int              eventqueueindexub;  /**< array position in event queue of upper bound change event, or -1 */
   int              parentvarssize;     /**< available slots in parentvars array */
   int              nparentvars;        /**< number of parent variables in aggregation tree (used slots of parentvars) */
   int              nuses;              /**< number of times, this variable is referenced */
   int              nlocksdown;         /**< number of locks for rounding down; if zero, rounding down is always feasible */
   int              nlocksup;           /**< number of locks for rounding up; if zero, rounding up is always feasible */
   unsigned int     inferdepth:16;      /**< depth in the tree, where this bound change took place */
   unsigned int     infernum:15;        /**< bound change index for each node representing the order of changes */
   unsigned int     removeable:1;       /**< TRUE iff var's column is removeable from the LP (due to aging or cleanup) */
   unsigned int     vartype:2;          /**< type of variable: binary, integer, implicit integer, continous */
   unsigned int     varstatus:3;        /**< status of variable: original, transformed, column, fixed, aggregated */
};




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
   TREE*            tree,               /**< branch-and-bound tree */
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
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** frees domain change data */
extern
RETCODE SCIPdomchgFree(
   DOMCHG**         domchg,             /**< pointer to domain change */
   MEMHDR*          memhdr              /**< block memory */
   );

/** converts a dynamic domain change data into a static one, using less memory than for a dynamic one */
extern
RETCODE SCIPdomchgMakeStatic(
   DOMCHG**         domchg,             /**< pointer to domain change data */
   MEMHDR*          memhdr              /**< block memory */
   );

/** applies domain change */
extern
RETCODE SCIPdomchgApply(
   DOMCHG*          domchg,             /**< domain change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
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
   TREE*            tree,               /**< branch-and-bound tree */
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
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
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
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
   Bool             removeable          /**< is var's column removeable from the LP (due to aging or cleanup)? */
   );

/** creates and captures a loose variable belonging to the transformed problem */
extern
RETCODE SCIPvarCreateTransformed(
   VAR**            var,                /**< pointer to variable data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   const char*      name,               /**< name of variable */
   Real             lb,                 /**< lower bound of variable */
   Real             ub,                 /**< upper bound of variable */
   Real             obj,                /**< objective function value */
   VARTYPE          vartype,            /**< type of variable */
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

/** increases lock number for rounding down; tells variable, that rounding its value down will make the solution
 *  infeasible
 */
extern
void SCIPvarForbidRoundDown(
   VAR*             var                 /**< problem variable */
   );

/** increases lock number for rounding up; tells variable, that rounding its value up will make the solution infeasible */
extern
void SCIPvarForbidRoundUp(
   VAR*             var                 /**< problem variable */
   );

/** increases lock number for rounding down and up; tells variable, that rounding value in either direction will make
 *  the solution infeasible
 */
extern
void SCIPvarForbidRound(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding down; cancels a prior forbidRoundDown() */
extern
void SCIPvarAllowRoundDown(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding up; cancels a prior forbidRoundUp() */
extern
void SCIPvarAllowRoundUp(
   VAR*             var                 /**< problem variable */
   );

/** decreases lock number for rounding down & up; cancels a prior forbidRound() */
extern
void SCIPvarAllowRound(
   VAR*             var                 /**< problem variable */
   );

/** gets number of locks for rounding down */
extern
int SCIPvarGetNLocksDown(
   VAR*             var                 /**< problem variable */
   );

/** gets number of locks for rounding up */
extern
int SCIPvarGetNLocksUp(
   VAR*             var                 /**< problem variable */
   );

/** is it possible, to round variable down and stay feasible? */
extern
Bool SCIPvarMayRoundDown(
   VAR*             var                 /**< problem variable */
   );

/** is it possible, to round variable up and stay feasible? */
extern
Bool SCIPvarMayRoundUp(
   VAR*             var                 /**< problem variable */
   );

/** copies original variable into loose transformed variable, that is captured */
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
   STAT*            stat                /**< problem statistics */
   );

/** converts variable into fixed variable */
extern
RETCODE SCIPvarFix(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             fixedval,           /**< value to fix variable at */
   Bool*            infeasible          /**< pointer to store whether the fixing is infeasible */
   );

/** converts variable into aggregated variable */
extern
RETCODE SCIPvarAggregate(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             aggvar,             /**< variable $y$ in aggregation $x = a*y + c$ */
   Real             scalar,             /**< multiplier $a$ in aggregation $x = a*y + c$ */
   Real             constant,           /**< constant shift $c$ in aggregation $x = a*y + c$ */
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
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   int              naggvars,           /**< number $n$ of variables in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   VAR**            aggvars,            /**< variables $y_i$ in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   Real*            scalars,            /**< multipliers $a_i$ in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   Real             constant,           /**< constant shift $c$ in aggregation $x = a_1*y_1 + ... + a_n*y_n + c$ */
   Bool*            infeasible          /**< pointer to store whether the aggregation is infeasible */
   );

/** gets negated variable x' = offset - x of problem variable x, where offset is fixed to lb + ub when the negated
 *  variable is created
 */
extern
RETCODE SCIPvarGetNegated(
   VAR*             var,                /**< problem variable to negate */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR**            negvar              /**< pointer to store the negated variable */
   );

/** gets the negation variable x of a negated variable x' = offset - x */
extern
VAR* SCIPvarGetNegationVar(
   VAR*             var                 /**< negated problem variable */
   );

/** gets the negation offset of a negated variable x' = offset - x */
extern
Real SCIPvarGetNegationConstant(
   VAR*             var                 /**< negated problem variable */
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
   TREE*            tree,               /**< branch-and-bound tree */
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
   TREE*            tree,               /**< branch-and-bound tree */
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
   TREE*            tree,               /**< branch-and-bound tree */
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
   TREE*            tree,               /**< branch-and-bound tree */
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
   TREE*            tree,               /**< branch-and-bound tree */
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

/** compares the index of two variables, returns -1 if first is smaller than, and +1 if first is greater than second
 *  variable index; returns 0 if both indices are equal, which means both variables are equal
 */
extern
int SCIPvarCmp(
   VAR*             var1,               /**< first problem variable */
   VAR*             var2                /**< second problem variable */
   );

/** gets corresponding active problem variable of a variable */
extern
VAR* SCIPvarGetProbvar(
   VAR*             var                 /**< problem variable */
   );

/** transforms given variable, boundtype and bound to the corresponding active variable values */
extern
RETCODE SCIPvarGetProbvarBound(
   VAR**            var,                /**< pointer to problem variable */
   Real*            bound,              /**< pointer to bound value to transform */
   BOUNDTYPE*       boundtype           /**< pointer to type of bound: lower or upper bound */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** get name of variable */
extern
const char* SCIPvarGetName(
   VAR*             var                 /**< problem variable */
   );

/** gets status of variable */
extern
VARSTATUS SCIPvarGetStatus(
   VAR*             var                 /**< problem variable */
   );

/** returns whether the variable belongs to the original problem */
extern
Bool SCIPvarIsOriginal(
   VAR*             var                 /**< problem variable */
   );

/** returns whether the variable belongs to the transformed problem */
extern
Bool SCIPvarIsTransformed(
   VAR*             var                 /**< problem variable */
   );

/** returns whether the variable was created by negation of a different variable */
extern
Bool SCIPvarIsNegated(
   VAR*             var                 /**< problem variable */
   );

/** gets type of variable */
extern
VARTYPE SCIPvarGetType(
   VAR*             var                 /**< problem variable */
   );

/** gets unique index of variable */
extern
int SCIPvarGetIndex(
   VAR*             var                 /**< problem variable */
   );

/** gets position of variable in problem, or -1 if variable is not active */
extern
int SCIPvarGetProbIndex(
   VAR*             var                 /**< problem variable */
   );

/** gets column of COLUMN variable */
extern
COL* SCIPvarGetCol(
   VAR*             var                 /**< problem variable */
   );

/** gets aggregation variable $y$ of an aggregated variable $x = a*y + c$ */
extern
VAR* SCIPvarGetAggrVar(
   VAR*             var                 /**< problem variable */
   );

/** gets aggregation scalar $a$ of an aggregated variable $x = a*y + c$ */
extern
Real SCIPvarGetAggrScalar(
   VAR*             var                 /**< problem variable */
   );

/** gets aggregation constant $c$ of an aggregated variable $x = a*y + c$ */
extern
Real SCIPvarGetAggrConstant(
   VAR*             var                 /**< problem variable */
   );

/** gets objective function value of variable */
extern
Real SCIPvarGetObj(
   VAR*             var                 /**< problem variable */
   );

/** gets global lower bound of variable */
extern
Real SCIPvarGetLbGlobal(
   VAR*             var                 /**< problem variable */
   );

/** gets global upper bound of variable */
extern
Real SCIPvarGetUbGlobal(
   VAR*             var                 /**< problem variable */
   );

/** gets current lower bound of variable */
extern
Real SCIPvarGetLbLocal(
   VAR*             var                 /**< problem variable */
   );

/** gets current upper bound of variable */
extern
Real SCIPvarGetUbLocal(
   VAR*             var                 /**< problem variable */
   );

/** gets inference constraint of variable (constraint that deduced the current assignment), or NULL */
extern
CONS* SCIPvarGetInferCons(
   VAR*             var                 /**< problem variable */
   );

/** gets inference variable of variable (variable that was assigned: parent of var, or var itself), or NULL */
extern
VAR* SCIPvarGetInferVar(
   VAR*             var                 /**< problem variable */
   );

/** gets inference depth level of variable */
extern
int SCIPvarGetInferDepth(
   VAR*             var                 /**< problem variable */
   );

/** gets inference number of variable (inference index in variable's inference depth level) */
extern
int SCIPvarGetInferNum(
   VAR*             var                 /**< problem variable */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPvarGetName(var)             (var)->name
#define SCIPvarGetStatus(var)           (VARSTATUS)((var)->varstatus)
#define SCIPvarIsOriginal(var)          ((var)->varstatus == SCIP_VARSTATUS_ORIGINAL \
      || ((var)->varstatus == SCIP_VARSTATUS_NEGATED && (var)->negatedvar->varstatus == SCIP_VARSTATUS_ORIGINAL))
#define SCIPvarIsTransformed(var)       ((var)->varstatus != SCIP_VARSTATUS_ORIGINAL \
      && ((var)->varstatus != SCIP_VARSTATUS_NEGATED || (var)->negatedvar->varstatus != SCIP_VARSTATUS_ORIGINAL))
#define SCIPvarIsNegated(var)           ((var)->varstatus == SCIP_VARSTATUS_NEGATED)
#define SCIPvarGetType(var)             ((VARTYPE)((var)->vartype))
#define SCIPvarGetIndex(var)            (var)->index
#define SCIPvarGetProbIndex(var)        (var)->probindex
#define SCIPvarGetCol(var)              (var)->data.col
#define SCIPvarGetAggrVar(var)          (var)->data.aggregate.var
#define SCIPvarGetAggrScalar(var)       (var)->data.aggregate.scalar
#define SCIPvarGetAggrConstant(var)     (var)->data.aggregate.constant
#define SCIPvarGetObj(var)              (var)->obj
#define SCIPvarGetLbGlobal(var)         (var)->glbdom.lb
#define SCIPvarGetUbGlobal(var)         (var)->glbdom.ub
#define SCIPvarGetLbLocal(var)          (var)->actdom.lb
#define SCIPvarGetUbLocal(var)          (var)->actdom.ub
#define SCIPvarGetInferCons(var)        (var)->infercons
#define SCIPvarGetInferVar(var)         (var)->infervar
#define SCIPvarGetInferDepth(var)       ((int)((var)->inferdepth))
#define SCIPvarGetInferNum(var)         ((int)((var)->infernum))

#endif

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

/** gets best local bound of variable with respect to the objective function */
extern
Real SCIPvarGetBestBound(
   VAR*             var                 /**< problem variable */
   );

/** gets type (lower or upper) of best bound of variable with respect to the objective function */
extern
BOUNDTYPE SCIPvarGetBestBoundType(
   VAR*             var                 /**< problem variable */
   );

/** gets primal LP solution value of variable */
extern
Real SCIPvarGetLPSol(
   VAR*             var                 /**< problem variable */
   );

/** gets pseudo solution value of variable at actual node */
extern
Real SCIPvarGetPseudoSol(
   VAR*             var                 /**< problem variable */
   );

/** gets solution value of variable at actual node: if LP was solved at the node, the method returns the LP primal
 *  solution value, otherwise the pseudo solution
 */
extern
Real SCIPvarGetSol(
   VAR*             var,                /**< problem variable */
   TREE*            tree                /**< branch-and-bound tree */
   );

/** resolves variable to columns and adds them with the coefficient to the row */
extern
RETCODE SCIPvarAddToRow(
   VAR*             var,                /**< problem variable */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
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
